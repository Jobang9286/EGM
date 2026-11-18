#!/usr/bin/env bash
PY="$HOME/.conda/envs/celloracle_env/bin/python"
SCRIPT="$(dirname "$0")/Fig3N_SuppFig9EF_plot_ko_flows_v2.py"

# -> Fig 3N Srebf2 KO, Nr4a1 KO (+ Control); Supp 9E Nfe2l1 KO
$PY $SCRIPT Srebf2 Srebf1 Nr4a1 Foxo1 Tfe3 Nr2f2 Nfe2l1 --tag lipid7 --scale-ref pool

# (same TFs at another pooled scale)
$PY $SCRIPT Srebf2 Nr4a1 Nfe2l2 Tfe3 Atf4 E2f1 Klf5 Atf6 --tag set8 --scale-ref pool

# -> Supp 9E Tfap2b KO, Klf5 KO
$PY $SCRIPT Srebf2 Nr4a1 Nfe2l2 Tfe3 Atf4 E2f1 Klf5 Atf6 Tfap2b Foxo1 Nr2f2 --tag set11 --scale-ref pool

# -> Fig 3N Mitf KO, Egr2 KO; Supp 9F Fli1 KO
$PY $SCRIPT Srebf2 Mitf Egr2 Fli1 Pbx1 Mbd2 Mafk --tag ctrl --scale-ref pool

