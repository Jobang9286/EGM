from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from mcp.server.fastmcp import FastMCP


def _import_rsem(
    samples: pd.DataFrame,
    rsem_dir: Path,
    suffix: str = ".genes.results",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Combine RSEM .genes.results files into gene × sample matrices."""
    counts_list, tpm_list = [], []
    for sample_id in samples["sample"]:
        path = rsem_dir / f"{sample_id}{suffix}"
        if not path.exists():
            raise FileNotFoundError(f"RSEM output not found: {path}")
        df = pd.read_csv(path, sep="\t")
        counts_list.append(df.set_index("gene_id")["expected_count"].rename(sample_id))
        tpm_list.append(df.set_index("gene_id")["TPM"].rename(sample_id))
    return pd.concat(counts_list, axis=1), pd.concat(tpm_list, axis=1)


def register_pipeline_tools(mcp: FastMCP, executor):
    q = executor.quote

    # ── Step 1: trim_galore + fastqc ─────────────────────────────
    @mcp.tool()
    def step1_preprocessing(
        project_dir: str,
        sample: str,
        num_threads: int = 0,
        trim_galore_bin: str = "trim_galore",
    ) -> str:
        """Step 1: Trim adapters with trim_galore and run FastQC.

        Auto-detects PE/SE and naming convention (`_1/_2`, `_R1/_R2`,
        `_R1_001/_R2_001`, bare `.fastq.gz`). Non-canonical inputs are
        symlinked to `{sample}_1/_2.fastq.gz` (PE) or `{sample}.fastq.gz`
        (SE) so downstream steps see canonical names.

        Args:
            project_dir:    Absolute path containing raw FASTQs.
            sample:         Sample prefix.
            num_threads:    CPU threads. 0 = auto-detect.
            trim_galore_bin: Path to trim_galore executable.
        """
        safe_dir = str(executor._resolve_project_dir(project_dir))
        if executor.is_done(safe_dir, sample, "step1"):
            return f"[SKIP] step1 already completed for {sample}"
        try:
            r1, r2, paired = _find_fastqs(safe_dir, sample)
        except FileNotFoundError as e:
            return f"[ERROR] {e}"
        _canonicalize_inputs(safe_dir, sample, r1, r2, paired)

        res = executor.detect_resources()
        t = num_threads or res["cores"]
        s = q(sample); d = q(safe_dir); tg = q(trim_galore_bin)
        mode = "paired" if paired else "single"
        header = f"""
set -euo pipefail
cd {d}
mkdir -p trim_galore fastqc
echo '[step1] mode={mode} sample={sample}'
"""
        if paired:
            body = f"""
{tg} --paired --gzip -o trim_galore/ {s}_1.fastq.gz {s}_2.fastq.gz
fastqc -t {t} --nogroup -o fastqc/ trim_galore/{s}_1_val_1.fq.gz
fastqc -t {t} --nogroup -o fastqc/ trim_galore/{s}_2_val_2.fq.gz
"""
        else:
            body = f"""
{tg} --gzip -o trim_galore/ {s}.fastq.gz
fastqc -t {t} --nogroup -o fastqc/ trim_galore/{s}_trimmed.fq.gz
"""
        out = executor.run(header + body, cwd=safe_dir)
        if executor.is_success(out):
            executor.mark_done(safe_dir, sample, "step1")
        return out

    # ── Step 2: STAR alignment ───────────────────────────────────
    @mcp.tool()
    def step2_star(
        project_dir: str,
        sample: str,
        ref_dir: str = "",
        num_threads: int = 0,
        star_bin: str = "STAR",
    ) -> str:
        """Step 2: STAR alignment → TranscriptomeSAM.

        Args:
            project_dir: Project directory with trim_galore/ outputs.
            sample:      Sample prefix.
            ref_dir:     STAR genomeDir. Default: Reference/RNA_hg38/.
            num_threads: CPU threads. 0 = auto-detect.
            star_bin:    Path to STAR executable.
        """
        safe_dir = str(executor._resolve_project_dir(project_dir))
        if executor.is_done(safe_dir, sample, "step2"):
            return f"[SKIP] step2 already completed for {sample}"
        try:
            tr1, tr2, paired = _find_trimmed(safe_dir, sample)
        except FileNotFoundError as e:
            return f"[ERROR] {e}"

        res = executor.detect_resources()
        t = num_threads or res["cores"]
        rd = ref_dir or "Reference/RNA_hg38/"
        s = q(sample); d = q(safe_dir)
        # Relative path from project_dir (trim_galore/XYZ) for cleaner command
        rel_r1 = "trim_galore/" + Path(tr1).name
        rel_r2 = "trim_galore/" + Path(tr2).name if paired else None
        reads_in = f"{rel_r1} {rel_r2}" if paired else rel_r1
        cmds = f"""
set -euo pipefail
cd {d}
echo '[step2] mode={"paired" if paired else "single"} sample={sample}'
# Resilience: STAR refuses to run if outTmpDir exists. Remove stale temp + partial outputs.
rm -rf TEMP_{s}
rm -f {s}Aligned.out.sam {s}Aligned.toTranscriptome.out.bam \
      {s}Log.final.out {s}Log.out {s}Log.progress.out {s}SJ.out.tab
{q(star_bin)} --runThreadN {t} --genomeDir {q(rd)} \
  --readFilesIn {reads_in} \
  --readFilesCommand zcat --outTmpDir TEMP_{s} \
  --outFileNamePrefix {s} --quantMode TranscriptomeSAM
rm -f {s}Aligned.out.sam
"""
        out = executor.run(cmds, cwd=safe_dir)
        if not executor.is_success(out):
            return out

        # Sanity: STAR must produce a complete transcriptome BAM and final log
        bam_out = f"{safe_dir}/{sample}Aligned.toTranscriptome.out.bam"
        star_log = f"{safe_dir}/{sample}Log.final.out"
        if not Path(star_log).exists():
            return out + f"\n[ERROR] STAR finished but {star_log} missing. Not marking step2 done."
        bam_err = _bam_quickcheck(executor, bam_out)
        if bam_err:
            return out + f"\n[ERROR] {bam_err}. Not marking step2 done."

        executor.mark_done(safe_dir, sample, "step2")
        return out

    # ── Step 3: RSEM quantification ──────────────────────────────
    @mcp.tool()
    def step3_rsem(
        project_dir: str,
        sample: str,
        ref_dir: str = "",
        genome: str = "GRCh38",
        num_threads: int = 0,
        rsem_bin: str = "rsem-calculate-expression",
    ) -> str:
        """Step 3: RSEM quantification → {sample}.genes.results + isoforms.results.

        Args:
            project_dir: Project directory.
            sample:      Sample prefix.
            ref_dir:     Reference directory.
            genome:      RSEM reference prefix (e.g. 'GRCh38').
            num_threads: CPU threads. 0 = auto-detect.
            rsem_bin:    Path to rsem-calculate-expression.
        """
        safe_dir = str(executor._resolve_project_dir(project_dir))
        if executor.is_done(safe_dir, sample, "step3"):
            return f"[SKIP] step3 already completed for {sample}"
        try:
            _, _, paired = _find_trimmed(safe_dir, sample)
        except FileNotFoundError as e:
            return f"[ERROR] {e}"

        # Validate STAR output is complete before running RSEM
        bam_in = f"{safe_dir}/{sample}Aligned.toTranscriptome.out.bam"
        star_log = f"{safe_dir}/{sample}Log.final.out"
        if not Path(star_log).exists():
            return f"[ERROR] STAR Log.final.out missing for '{sample}' — step2 incomplete or failed."
        bam_err = _bam_quickcheck(executor, bam_in)
        if bam_err:
            return f"[ERROR] {bam_err}. Re-run step2 first."

        res = executor.detect_resources()
        t = num_threads or res["cores"]
        rd = ref_dir or "Reference/RNA_hg38/"
        ref_prefix = str(Path(rd) / genome)  # correct path join regardless of trailing slash
        s = q(sample); d = q(safe_dir)
        pe_flag = "--paired-end" if paired else ""
        cmds = f"""
set -euo pipefail
cd {d}
echo '[step3] mode={"paired" if paired else "single"} sample={sample}'
# Resilience: clear stale RSEM intermediates from interrupted prior run
rm -rf {s}.temp {s}.stat
rm -f {s}.genes.results {s}.isoforms.results {s}.transcript.bam
{q(rsem_bin)} {pe_flag} --bam -p {t} \
  {s}Aligned.toTranscriptome.out.bam {q(ref_prefix)} {s}
"""
        out = executor.run(cmds, cwd=safe_dir)
        if not executor.is_success(out):
            return out

        # Output sanity: RSEM with empty/corrupt input still exits 0; reject empty outputs
        genes_file = f"{safe_dir}/{sample}.genes.results"
        total = _rsem_genes_total_counts(genes_file)
        if total < 1000:
            return (out + f"\n[ERROR] RSEM produced near-empty output "
                          f"(total expected_count={total:.0f}). Not marking step3 done. "
                          f"Likely the input BAM was truncated.")

        executor.mark_done(safe_dir, sample, "step3")
        return out + f"\n[step3] OK: {sample} total expected_count={total:.0f}"

    # ── Step 4: Load RSEM counts (tximport equivalent) ───────────
    @mcp.tool()
    def step4_load_counts(
        project_dir: str,
        mapfile: str,
        rsem_dir: str = "",
        suffix: str = ".genes.results",
    ) -> str:
        """Step 4: Combine per-sample RSEM outputs into counts/TPM matrices.

        Args:
            project_dir: Project directory (outputs written here).
            mapfile:     TSV path with columns 'sample' and 'TREATMENT'.
            rsem_dir:    Directory with .genes.results files. Defaults to project_dir.
            suffix:      RSEM filename suffix.
        """
        pd_dir = executor._resolve_project_dir(project_dir)
        mapfile_path = _resolve_input_path(executor, mapfile)
        samples = pd.read_csv(mapfile_path, sep="\t")
        if not {"sample", "TREATMENT"} <= set(samples.columns):
            return "[ERROR] mapfile must have columns 'sample' and 'TREATMENT'"
        samples["TREATMENT"] = samples["TREATMENT"].str.lower()

        rd = _resolve_input_path(executor, rsem_dir) or pd_dir
        counts, tpm = _import_rsem(samples, rd, suffix)
        counts.to_csv(pd_dir / "txi.rsem.counts.txt", sep="\t")
        tpm.to_csv(pd_dir / "txi.rsem.abundance.TPM.txt", sep="\t")

        executor.mark_done(str(pd_dir), "_global", "step4")
        return json.dumps({
            "counts": str(pd_dir / "txi.rsem.counts.txt"),
            "tpm":    str(pd_dir / "txi.rsem.abundance.TPM.txt"),
            "shape":  {"genes": counts.shape[0], "samples": counts.shape[1]},
            "samples": list(counts.columns),
        }, indent=2)

