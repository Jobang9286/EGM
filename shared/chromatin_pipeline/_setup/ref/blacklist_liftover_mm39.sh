#!/bin/bash
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"

cd $ROOT/4FP_blood/_setup/ref

## ---- decompress the mm10 list --------------------------------------------------
zcat mm10-blacklist.v2.bed.gz > mm10-blacklist.v2.bed

## ---- keep chrom/start/end, lift to mm39, sort ----------------------------------
cut -f1-3 mm10-blacklist.v2.bed | sort -k1,1 -k2,2n > mm10-bl3.bed
./liftOver mm10-bl3.bed mm10ToMm39.over.chain.gz mm39-blacklist.v2.bed mm39-blacklist.unmapped 2>&1 | tail -3
echo "mm10 in : $(wc -l < mm10-bl3.bed)"
echo "mm39 out: $(wc -l < mm39-blacklist.v2.bed)"
echo "unmapped: $(grep -vc '^#' mm39-blacklist.unmapped 2>/dev/null || echo 0)"
sort -k1,1 -k2,2n mm39-blacklist.v2.bed -o mm39-blacklist.v2.bed
