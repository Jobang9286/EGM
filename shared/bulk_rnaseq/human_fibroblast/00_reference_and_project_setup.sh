#!/bin/bash
ROOT=$PWD
REF_SRC=$ROOT/GRCh38

## ---- A. project layout, sample sheet, tool shims, GTF -------------------------
set -e
BASE=$ROOT/4_human_egm
mkdir -p $BASE/project $BASE/_bin $BASE/_scripts $BASE/_logs $BASE/project/Reference/hg38
cd $BASE/project
for s in Control-1 Control-2 Control-3 EGM-1 EGM-2 EGM-3; do
  ln -sfn ../Raw_data/$s/${s}_1.fastq.gz ${s}_1.fastq.gz
  ln -sfn ../Raw_data/$s/${s}_2.fastq.gz ${s}_2.fastq.gz
done
printf 'sample\tTREATMENT\nControl-1\tuntreated\nControl-2\tuntreated\nControl-3\tuntreated\nEGM-1\ttreated\nEGM-2\ttreated\nEGM-3\ttreated\n' > mapfile.tsv
# tool shim dir: trim_galore wrapper (multi-core + explicit cutadapt), cutadapt from chipseq env
cat > $BASE/_bin/trim_galore <<EOF
#!/bin/bash
exec $HOME/.conda/envs/rnaseq/bin/trim_galore --cores 4 --path_to_cutadapt $HOME/.conda/envs/chipseq/bin/cutadapt "\$@"
EOF
chmod +x $BASE/_bin/trim_galore
ln -sfn $HOME/.conda/envs/chipseq/bin/cutadapt $BASE/_bin/cutadapt
# GTF inside the sandbox (hard link, same filesystem) for annotate_symbols
ln -f $REF_SRC/Homo_sapiens.GRCh38.110.gtf $BASE/project/Reference/hg38/Homo_sapiens.GRCh38.110.gtf 2>/dev/null || cp $REF_SRC/Homo_sapiens.GRCh38.110.gtf $BASE/project/Reference/hg38/

## ---- B. RSEM reference -----------------------------------------------------------
mkdir -p $ROOT/human_ref/RNA_hg38 $ROOT/human_ref/logs
cd $ROOT/human_ref
export PATH=$HOME/.conda/envs/rnaseq/bin:$PATH
setsid nohup rsem-prepare-reference --num-threads 16 --gtf $REF_SRC/Homo_sapiens.GRCh38.110.gtf $REF_SRC/Homo_sapiens.GRCh38.dna.primary_assembly.fa $ROOT/human_ref/RNA_hg38/GRCh38 > logs/rsem_build.log 2>&1 < /dev/null &

## ---- C. STAR index ---------------------------------------------------------------
cd $ROOT/human_ref && export PATH=$HOME/.conda/envs/rnaseq/bin:$PATH && setsid nohup STAR --runThreadN 32 --runMode genomeGenerate --genomeDir $ROOT/human_ref/RNA_hg38 --genomeFastaFiles $REF_SRC/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile $REF_SRC/Homo_sapiens.GRCh38.110.gtf --sjdbOverhang 100 --limitGenomeGenerateRAM 128000000000 --outFileNamePrefix $ROOT/human_ref/logs/star_ > logs/star_build.log 2>&1 < /dev/null &

## ---- D. gene id -> symbol table --------------------------------------------------
# genes without gene_name keep their Ensembl id as the name (20,393 of 62,754)
cd $ROOT/4_human_egm/project/Reference/hg38 && { printf 'Gene.ID\tGene.Name\n'; awk -F'\t' '$3=="gene"{ if (match($9,/gene_id "[^"]+"/)) id=substr($9,RSTART+9,RLENGTH-10); else next; n=""; if (match($9,/gene_name "[^"]+"/)) n=substr($9,RSTART+11,RLENGTH-12); print id"\t"(n==""?id:n)}' Homo_sapiens.GRCh38.110.gtf; } > gene_annotation.tsv

## ---- E. launch trimming + FastQC (phase step1 of 01_run_rna_mcp.py) ------------
# (01_run_rna_mcp.py is installed as $BASE/_scripts/run_rna_mcp.py; 02_orchestrate.sh
#  was then launched to chain the remaining phases)
cd $ROOT && setsid nohup $HOME/RNA_mcp/.venv/bin/python -I $BASE/_scripts/run_rna_mcp.py step1 > $BASE/_logs/driver_step1.log 2>&1 < /dev/null &
