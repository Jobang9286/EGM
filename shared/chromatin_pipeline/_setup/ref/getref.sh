#!/bin/bash
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
set -x
cd $ROOT/4FP_blood/_setup/ref
curl -sL -o GRCm39.zip https://genome-idx.s3.amazonaws.com/bt/GRCm39.zip &
curl -sL -o mm10.zip    https://genome-idx.s3.amazonaws.com/bt/mm10.zip &
curl -sL -o mm10-blacklist.v2.bed.gz https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz
curl -sL -o mm10ToMm39.over.chain.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz
curl -sL -o liftOver https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver && chmod +x liftOver
wait
ls -lh
