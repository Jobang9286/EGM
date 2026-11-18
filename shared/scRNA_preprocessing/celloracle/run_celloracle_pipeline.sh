#!/usr/bin/env bash
ORACLE="EGM_multiflex/oracle"
PY="$HOME/.conda/envs/celloracle_env/bin/python"
CODE="$(cd "$(dirname "$0")" && pwd)"
cd "$ORACLE"

# ---- GRN (phase 1), centrality (phase 2) ----
$PY $CODE/01b_celloracle_phase1_max.py > phase1_max_$(date +%Y%m%d_%H%M%S).log 2>&1
$PY $CODE/02_celloracle_phase2_v2.py

# ---- 34 curated lipid TFs -> batch 1 (cosine scoring) ----
$PY $CODE/04a_make_lipid_tf_list.py
nohup $PY $CODE/05_celloracle_ko_batch.py --tag lipid --tf-file lipid_tf_list.txt --n-jobs 16 > batch_lipid.log 2>&1 &

# ---- GO-derived lipid regulators -> batches 2a / 2b ----
Rscript $CODE/04b_lipid_go_tfs.R
$PY $CODE/04c_make_lipid2_tf_list.py
head -25 lipid2_tf_list.txt > lipid2a_tf_list.txt; tail -n +26 lipid2_tf_list.txt > lipid2b_tf_list.txt
setsid nohup $PY $CODE/05_celloracle_ko_batch.py --tag lipid2a --no-ref --tf-file lipid2a_tf_list.txt --n-jobs 16 </dev/null >batch_lipid2a.log 2>&1 &
setsid nohup $PY $CODE/05_celloracle_ko_batch.py --tag lipid2b --no-ref --tf-file lipid2b_tf_list.txt --n-jobs 16 </dev/null >batch_lipid2b.log 2>&1 &
# lipid2b was OOM-killed; restarted alone
mv batch_lipid2b.log batch_lipid2b_attempt1.log
setsid nohup $PY $CODE/05_celloracle_ko_batch.py --tag lipid2b --no-ref --tf-file lipid2b_tf_list.txt --n-jobs 16 </dev/null >batch_lipid2b.log 2>&1 &
wait

# ---- batch-1 annotation tables ----
$PY $CODE/06_finalize_lipid_ko.py
$PY $CODE/07_finalize_lipid_ko_curated.py
$PY $CODE/08_build_supp_all84_table.py

# ---- CellOracle perturbation score for all 84 lipid-associated TFs ----
$PY $CODE/09a_make_v2_tf_lists.py
setsid nohup $PY $CODE/09b_celloracle_ko_batch_v2.py --tag v2a --no-ref --tf-file v2a_tf_list.txt --n-jobs 16 </dev/null >batch_v2a.log 2>&1 &
setsid nohup $PY $CODE/09b_celloracle_ko_batch_v2.py --tag v2b --no-ref --tf-file v2b_tf_list.txt --n-jobs 16 </dev/null >batch_v2b.log 2>&1 &
wait

# ---- BH correction and final lipid-TF table (input of Fig 3M) ----
$PY $CODE/10_finalize_ps_v2.py

# ---- non-lipid negative-control TFs (Fig 3M grey bars Mbd2 = negctrl2, Fli1 = negctrl3) ----
printf "Stat1\nStat2\nIrf7\nRela\nFoxm1\nE2f6\nTead1\nTead2\nCtcf\nEts1\nElf1\nTwist2\n" > negctrl_tf_list.txt
setsid nohup $PY $CODE/09b_celloracle_ko_batch_v2.py --tag negctrl --no-ref --tf-file negctrl_tf_list.txt --n-jobs 32 </dev/null >batch_negctrl.log 2>&1 &
printf "Mbd2\nSmarcc2\nSetdb1\nTrim28\nHdac2\nSmc3\n" > negctrl2_tf_list.txt
setsid nohup $PY $CODE/09b_celloracle_ko_batch_v2.py --tag negctrl2 --no-ref --tf-file negctrl2_tf_list.txt --n-jobs 32 </dev/null >batch_negctrl2.log 2>&1 &
printf "Ets2\nElk3\nEtv5\nGabpa\nFli1\n" > negctrl3_tf_list.txt
setsid nohup $PY $CODE/09b_celloracle_ko_batch_v2.py --tag negctrl3 --no-ref --tf-file negctrl3_tf_list.txt --n-jobs 32 </dev/null >batch_negctrl3.log 2>&1 &
# (negctrl4 = Irf1 Irf3 Irf9 Rel Nfkb2 was stopped without results)
wait
# PS_negctrl_all.csv (used for Fig 3M)
$PY $CODE/11a_build_PS_negctrl_all.py

# ---- negctrl6: its Mafk KO vectors enter the pooled arrow scale of the ctrl flow plot (Fig 3N, Supp 9F) ----
printf "Tcf12\nHic1\nMafk\nZfp281\nPrrx1\nSmad4\nIrf9\nBach1\n" > negctrl6_tf_list.txt
setsid nohup $PY $CODE/09b_celloracle_ko_batch_v2.py --tag negctrl6 --no-ref --tf-file negctrl6_tf_list.txt --n-jobs 32 </dev/null >batch_negctrl6.log 2>&1 &
wait

# ---- copy the two tables to the local machine for Fig 3M ----
#   phase3_figs/PS_TF_KO_all.csv    -> data/PS_TF_KO_all.csv
#   phase3_figs/PS_negctrl_all.csv  -> data/PS_negctrl_all.csv
# Vector-field panels (Fig 3N, Supp Fig 9E-F): ../../../Figure3/celloracle/Fig3N_SuppFig9EF_run_flow_plots.sh
