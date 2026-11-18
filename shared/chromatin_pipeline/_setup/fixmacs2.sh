#!/bin/bash
set -x
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate chipseq
# The bioconda MACS2 binary was linked against a glibc that still exported
# __log_finite; RHEL 9 (glibc 2.34) does not, so the .so fails to load.
# Recompiling only MACS2 against this glibc fixes it; deps stay as wheels.
pip install --no-cache-dir --force-reinstall --no-binary macs2 macs2==2.2.9.1
python -c "from MACS2.IO import ScoreTrack; print('IMPORT_OK')"
macs2 --version && echo MACS2_FIXED
