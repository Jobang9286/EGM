#!/usr/bin/env python3
"""Drive the RNA_mcp pipeline tools for 4_human_egm (Control x3 vs EGM x3, GRCh38).

Usage: run_rna_mcp.py <phase> [phase ...]   phases: step1 | align | de | qc
Calls the registered MCP tools (step1..step8, run_full_pipeline, AI QC tools) directly.
"""
import json, os, sys, time, threading
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.getcwd()

BASE = f"{ROOT}/4_human_egm"
PROJ = f"{BASE}/project"
REF = f"{ROOT}/human_ref/RNA_hg38/"
GENOME = "GRCh38"
SAMPLES = ["Control-1", "Control-2", "Control-3", "EGM-1", "EGM-2", "EGM-3"]
MAPFILE = f"{PROJ}/mapfile.tsv"
GTF = f"{PROJ}/Reference/hg38/Homo_sapiens.GRCh38.110.gtf"
ANNOT = f"{PROJ}/Reference/hg38/gene_annotation.tsv"

os.environ["RNA_WORKDIR"] = BASE
os.environ["RNA_MCP_PYTHON"] = os.path.expanduser("~/RNA_mcp/.venv/bin/python")
os.environ["PATH"] = f"{BASE}/_bin:{os.path.expanduser('~/.conda/envs/rnaseq/bin')}:" + os.environ.get("PATH", "")

from rna_mcp.server import mcp  # noqa: E402  (executor reads RNA_WORKDIR at import)
TM = mcp._tool_manager
_lock = threading.Lock()


def log(msg):
    with _lock:
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def call(name, **kw):
    t0 = time.time()
    try:
        out = TM.get_tool(name).fn(**kw)
    except Exception as e:  # keep the driver alive, record the failure
        out = f"[EXCEPTION] {type(e).__name__}: {e}"
    log(f"{name} {json.dumps(kw, ensure_ascii=False)} -> {time.time()-t0:.0f}s\n{out}\n")
    return out


def phase_step1():
    with ThreadPoolExecutor(max_workers=len(SAMPLES)) as pool:
        futs = [pool.submit(call, "step1_preprocessing", project_dir=PROJ, sample=s, num_threads=8)
                for s in SAMPLES]
        for f in as_completed(futs):
            f.result()


def phase_align():
    call("run_full_pipeline", project_dir=PROJ, ref_dir=REF, genome=GENOME)


def phase_de():
    call("step4_load_counts", project_dir=PROJ, mapfile=MAPFILE)
    call("step5_deseq2", project_dir=PROJ, mapfile=MAPFILE, engine="r")
    call("step6_plot_ma", project_dir=PROJ)
    call("step7_plot_volcano", project_dir=PROJ)
    call("step8_plot_heatmap", project_dir=PROJ, mapfile=MAPFILE, annotation_file=ANNOT)
    call("annotate_symbols", project_dir=PROJ, gtf_file=GTF)


def phase_qc():
    for s in SAMPLES:
        call("check_qc", project_dir=PROJ, sample=s)
        call("detect_strandedness", project_dir=PROJ, sample=s)
    call("rna_status", project_dir=PROJ)


def phase_sample(sample):
    """step2 (STAR) + step3 (RSEM) for one sample, outside run_full_pipeline's lock."""
    call("step2_star", project_dir=PROJ, sample=sample, ref_dir=REF, num_threads=32)
    call("step3_rsem", project_dir=PROJ, sample=sample, ref_dir=REF, genome=GENOME, num_threads=32)


def phase_rsem_reverse():
    """step3 (RSEM) for all samples with --strandedness reverse via a wrapper rsem_bin."""
    rsem_bin = f"{BASE}/_bin/rsem-calculate-expression-reverse"
    with ThreadPoolExecutor(max_workers=len(SAMPLES)) as pool:
        futs = [pool.submit(call, "step3_rsem", project_dir=PROJ, sample=s, ref_dir=REF,
                            genome=GENOME, num_threads=20, rsem_bin=rsem_bin) for s in SAMPLES]
        for f in as_completed(futs):
            f.result()


PHASES = {"step1": phase_step1, "align": phase_align, "de": phase_de, "qc": phase_qc,
          "rsem_reverse": phase_rsem_reverse}

if __name__ == "__main__":
    if len(sys.argv) < 2 or any(p not in PHASES and not p.startswith("sample:") for p in sys.argv[1:]):
        sys.exit(f"usage: {sys.argv[0]} {{{'|'.join(PHASES)}|sample:<name>}} ...")
    log(f"tools available: {[t.name for t in TM.list_tools()]}")
    for ph in sys.argv[1:]:
        log(f"=== phase {ph} START ===")
        if ph.startswith("sample:"):
            phase_sample(ph.split(":", 1)[1])
        else:
            PHASES[ph]()
        log(f"=== phase {ph} END ===")
