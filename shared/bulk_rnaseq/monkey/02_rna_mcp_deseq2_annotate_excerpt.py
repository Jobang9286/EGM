"""Pipeline tools — RNA-seq shell steps + DE analysis + visualisation."""
from __future__ import annotations

import fcntl
import glob as g
import json
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from mcp.server.fastmcp import FastMCP
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats


def _resolve_input_path(executor, path: str) -> Path | None:
    """Resolve optional user-supplied file/dir paths under the configured work_dir."""
    if not path:
        return None
    return executor._resolve_project_dir(path)


import subprocess as _subproc
import tempfile as _tempfile
import shutil as _shutil


def _run_rscript(r_code: str, cwd: str | None = None) -> tuple[int, str, str]:
    """Run R code via Rscript. Returns (returncode, stdout, stderr)."""
    rscript = _shutil.which("Rscript") or "/usr/bin/Rscript"
    with _tempfile.NamedTemporaryFile("w", suffix=".R", delete=False) as fh:
        fh.write(r_code)
        path = fh.name
    try:
        proc = _subproc.run([rscript, path], cwd=cwd, capture_output=True, text=True)
        return proc.returncode, proc.stdout, proc.stderr
    finally:
        Path(path).unlink(missing_ok=True)


def _r_str(s: str) -> str:
    """Quote a string for safe embedding in R code."""
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def _run_deseq2_r(
    counts_file: Path,
    mapfile: Path,
    out_results: Path,
    out_norm: Path,
    reference: str = "untreated",
) -> dict:
    """Run R DESeq2 with default settings (independent filtering ON).

    Output format matches pydeseq2 columns: baseMean, log2FoldChange, lfcSE,
    stat, pvalue, padj — so downstream step6/7/8 plots work unchanged.
    """
    r_code = f"""
suppressPackageStartupMessages({{ library(DESeq2) }})
counts <- read.table({_r_str(str(counts_file))}, sep="\\t",
                    header=TRUE, row.names=1, check.names=FALSE)
samples <- read.table({_r_str(str(mapfile))}, header=TRUE, sep="\\t",
                     stringsAsFactors=FALSE)
samples$TREATMENT <- tolower(samples$TREATMENT)

# Reorder counts columns to match mapfile, drop any unmatched
common <- intersect(samples$sample, colnames(counts))
samples <- samples[match(common, samples$sample), ]
counts  <- counts[, common, drop=FALSE]
counts  <- round(as.matrix(counts))
mode(counts) <- "integer"

# Reference level first
samples$TREATMENT <- factor(samples$TREATMENT,
                            levels=unique(c({_r_str(reference.lower())},
                                            sort(unique(samples$TREATMENT)))))
samples$TREATMENT <- droplevels(samples$TREATMENT)
coldata <- data.frame(row.names=samples$sample, TREATMENT=samples$TREATMENT)

dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata,
                              design=~TREATMENT)
dds <- DESeq(dds, quiet=TRUE)
res <- results(dds, contrast=c("TREATMENT", "treated", {_r_str(reference.lower())}))
# Independent filtering applied by default; reference R DESeq2 behavior
res <- res[order(res$padj, na.last=TRUE), ]

resdf <- as.data.frame(res)
resdf <- cbind(gene_id=rownames(resdf), resdf)
write.table(resdf, {_r_str(str(out_results))}, sep=",", quote=FALSE,
            row.names=FALSE, col.names=TRUE)

# Save normalized counts (genes x samples)
norm <- counts(dds, normalized=TRUE)
normdf <- cbind(gene_id=rownames(norm), as.data.frame(norm))
write.table(normdf, {_r_str(str(out_norm))}, sep="\\t", quote=FALSE,
            row.names=FALSE, col.names=TRUE)

# Print summary as JSON-ish lines for caller
total <- nrow(resdf)
tested <- sum(!is.na(resdf$padj))
sig <- sum(!is.na(resdf$padj) & resdf$padj < 0.10 &
           abs(resdf$log2FoldChange) >= 1)
cat(sprintf("TOTAL=%d\\nTESTED=%d\\nSIGNIFICANT=%d\\n", total, tested, sig))
"""
    rc, stdout, err = _run_rscript(r_code)
    if rc != 0:
        raise RuntimeError(f"R DESeq2 failed: {err}")
    info = {"total_genes": 0, "tested_genes": 0, "significant": 0}
    for line in stdout.splitlines():
        if line.startswith("TOTAL="):
            info["total_genes"] = int(line.split("=", 1)[1])
        elif line.startswith("TESTED="):
            info["tested_genes"] = int(line.split("=", 1)[1])
        elif line.startswith("SIGNIFICANT="):
            info["significant"] = int(line.split("=", 1)[1])
    return info


def register_pipeline_tools(mcp: FastMCP, executor):
    # ── Step 5: DESeq2 DE (R or pydeseq2) ────────────────────────
    @mcp.tool()
    def step5_deseq2(
        project_dir: str,
        mapfile: str,
        counts_file: str = "",
        reference: str = "untreated",
        num_threads: int = 8,
        engine: str = "r",
        quant_source: str = "counts",
    ) -> str:
        """Step 5: DESeq2 differential expression analysis.

        engine='r' (default): R DESeq2 (reference implementation,
                              independent filtering ON by default).
        engine='python': pydeseq2 (legacy / pure-Python pipeline).

        Args:
            project_dir: Project directory (outputs written here).
            mapfile:     TSV with columns 'sample', 'TREATMENT' (treated/untreated).
            counts_file: Path to counts matrix. Defaults to txi.rsem.counts.txt.
            reference:   Reference condition name.
            num_threads: Parallel workers (only used by pydeseq2 engine).
            engine:      'r' or 'python'.
            quant_source: 'counts' (default, mirrors 1.DESeq2.R via
                          txi.rsem.counts.txt) or 'tpm' (mirrors
                          1.DESeq2_TPM.R via txi.rsem.abundance.TPM.txt;
                          rounded to integer inside DESeq2).
        """
        pd_dir = executor._resolve_project_dir(project_dir)
        _default_matrix = ("txi.rsem.abundance.TPM.txt"
                           if quant_source == "tpm" else "txi.rsem.counts.txt")
        counts_path = _resolve_input_path(executor, counts_file) or pd_dir / _default_matrix
        mapfile_path = _resolve_input_path(executor, mapfile)
        results_out = pd_dir / "DESeq2_results.csv"
        norm_out = pd_dir / "DESeq2_normalized_counts.tsv"

        if engine == "r":
            try:
                info = _run_deseq2_r(counts_path, mapfile_path, results_out, norm_out,
                                     reference=reference)
            except RuntimeError as e:
                return f"[ERROR] {e}"
            executor.mark_done(str(pd_dir), "_global", "step5")
            return json.dumps({
                "engine":          "r",
                "results":         str(results_out),
                "normalized":      str(norm_out),
                "total_genes":     info["total_genes"],
                "tested_genes":    info["tested_genes"],
                "significant":     info["significant"],
                "significant_thr": "padj<0.10 & |log2FC|>=1",
            }, indent=2)

        # Python (pydeseq2) engine — legacy
        counts = pd.read_csv(counts_path, sep="\t", index_col=0)
        samples = pd.read_csv(mapfile_path, sep="\t")
        samples["TREATMENT"] = samples["TREATMENT"].str.lower()
        err = _validate_deseq_inputs(list(counts.columns), samples, reference.lower())
        if err:
            return f"[ERROR] {err}"

        results, dds = _run_deseq2(counts, samples, reference, num_threads)
        results.to_csv(results_out)

        norm = pd.DataFrame(dds.layers["normed_counts"],
                            index=dds.obs_names, columns=dds.var_names).T
        norm.to_csv(norm_out, sep="\t")

        sig = results[(results["padj"] < 0.10) & (results["log2FoldChange"].abs() >= 1)]
        executor.mark_done(str(pd_dir), "_global", "step5")
        return json.dumps({
            "engine":          "python",
            "results":         str(results_out),
            "normalized":      str(norm_out),
            "total_genes":     len(results),
            "significant":     len(sig),
            "significant_thr": "padj<0.10 & |log2FC|>=1",
        }, indent=2)


# ---- gene ID -> symbol annotation ----------------------------------------------
def _resolve_input_path(executor, path: str) -> Path | None:
    """Resolve optional user-supplied file paths under the configured work_dir."""
    if not path:
        return None
    return executor._resolve_project_dir(path)


def register_ai_tools(mcp: FastMCP, executor):
    """Register 11 AI inference tools (9 adapted from ATAC + 2 RNA-specific)."""
    q = executor.quote
    # ── 11. Gene ID → Symbol Annotation ─────────────────────────
    @mcp.tool()
    def annotate_symbols(
        project_dir: str,
        results_file: str = "",
        gtf_file: str = "",
        out_name: str = "DESeq2_results_annotated.csv",
    ) -> str:
        """Map Ensembl gene IDs in DE results to gene symbols using a GTF.
        Looks up 'gene_id' → 'gene_name' attributes in the GTF.

        Args:
            project_dir:  Project directory.
            results_file: DE results CSV. Default: project_dir/DESeq2_results.csv.
            gtf_file:     GTF path. Default tries {project_dir}/Reference/hg38/*.gtf.
            out_name:     Output filename (written into project_dir).
        """
        pd_dir = executor._resolve_project_dir(project_dir)
        results_path = _resolve_input_path(executor, results_file) or pd_dir / "DESeq2_results.csv"
        if not results_path.exists():
            return json.dumps({"error": f"Results not found: {results_path}"})

        gtf_path = _resolve_input_path(executor, gtf_file)
        if gtf_path is None:
            cand = list((pd_dir / "Reference" / "hg38").glob("*.gtf")) + \
                   list(pd_dir.glob("*.gtf"))
            if not cand:
                return json.dumps({
                    "error": "GTF not specified and none found under project_dir. "
                             "Provide gtf_file argument."
                })
            gtf_path = cand[0]

        # Parse GTF attributes column (9th, tab-separated).
        # Store both versioned (ENSG00000000003.15) and unversioned
        # (ENSG00000000003) keys so GENCODE-style IDs match regardless of
        # whether the DE results carry a version suffix or not.
        id_to_name: dict[str, str] = {}
        with open(gtf_path) as f:
            for line in f:
                if line.startswith("#") or "\tgene\t" not in line:
                    continue
                attr = line.rstrip().split("\t")[8]
                gid = _gtf_attr(attr, "gene_id")
                gname = _gtf_attr(attr, "gene_name") or gid
                if gid:
                    id_to_name[gid] = gname
                    base = gid.split(".")[0]
                    if base not in id_to_name:
                        id_to_name[base] = gname

        res = pd.read_csv(results_path, index_col=0)
        res.insert(0, "gene_symbol", [
            id_to_name.get(g) or id_to_name.get(g.split(".")[0], g)
            for g in res.index
        ])
        out_path = pd_dir / out_name
        res.to_csv(out_path)
        mapped = sum(1 for g in res.index if g in id_to_name)
        return json.dumps({
            "output": str(out_path),
            "total_genes": len(res),
            "mapped_to_symbol": mapped,
            "gtf_used": str(gtf_path),
        }, indent=2)


def _gtf_attr(attr: str, key: str) -> str | None:
    """Extract a value like gene_id "ENSG..."; from a GTF attributes string."""
    needle = f'{key} "'
    i = attr.find(needle)
    if i < 0:
        return None
    i += len(needle)
    j = attr.find('"', i)
    return attr[i:j] if j > i else None
