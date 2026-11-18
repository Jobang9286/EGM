import os
import glob
import subprocess

READS_DIR = "39_fibroblasts/2_RNA_seq_tissue/Kidney_final"
OUTPUT_DIR = "39_fibroblasts/2_RNA_seq_tissue/Kidney_final/bam"
GENOME_DIR = "STAR_ref/GRCm39"

os.makedirs(OUTPUT_DIR, exist_ok=True)

samples = sorted(set([os.path.basename(f).split('_1.fastq.gz')[0] for f in glob.glob(f"{READS_DIR}/*_1.fastq.gz")]))


for sample in samples:
    R1 = os.path.join(READS_DIR, f"{sample}_1.fastq.gz")
    R2 = os.path.join(READS_DIR, f"{sample}_2.fastq.gz")
    PREFIX = os.path.join(OUTPUT_DIR, f"{sample}_")

    print(f"Processing sample {sample} with R1: {R1} and R2: {R2}")

    if not os.path.isfile(R1):
        print(f"Error: File {R1} does not exist!")
        continue
    if not os.path.isfile(R2):
        print(f"Error: File {R2} does not exist!")
        continue

    command = [
        "STAR",
        "--genomeDir", GENOME_DIR,
        "--readFilesIn", R1, R2,
        "--readFilesCommand", "zcat",
        "--runThreadN", "48",
        "--outFileNamePrefix", PREFIX,
        "--outSAMtype", "BAM", "SortedByCoordinate"
    ]

    try:
        subprocess.run(command, check=True)
        print(f"Finished processing {sample}")
    except subprocess.CalledProcessError as e:
        print(f"Error running STAR for sample {sample}: {e}")
