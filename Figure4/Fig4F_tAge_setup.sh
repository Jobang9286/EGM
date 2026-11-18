#!/usr/bin/env bash
set -e
HERE="$(pwd)"

## 1) tAge R package — install from the bundled source tarball (no GitHub auth needed)
tar xzf vendor/tAge_source.tar.gz -C vendor            # -> vendor/tAge-main
R CMD INSTALL vendor/tAge-main

## 2) Python venv (reticulate backend for tAge's .pkl ElasticNet models)
python3 -m venv .venv
.venv/bin/pip install --quiet --upgrade pip
.venv/bin/pip install --quiet -r vendor/requirements.txt   # joblib pandas scikit-learn numpy

## 3) Clock models: already provided in clocks/ (pre-downloaded from Zenodo 18763485).
##    To (re)download the full set instead:
##    Rscript -e 'Sys.setenv(RETICULATE_PYTHON=".venv/bin/python"); library(tAge);
##      cl<-list_clocks();
##      download_clocks(cl[grepl("EN_(Mortality|NormalizedAge|Chronoage)_Rodents",cl$filename),], "clocks")'

echo ""
echo "setup complete. Reproduce Fig 4F with:"
echo "  Rscript Fig4F_mortality_tAge_predict.R   # -> outputs/tage_mortality_persample.csv"
echo "  Rscript Fig4F_mortality_tAge_plot.R      # -> outputs/mouse_mortality_fig4C.png"
