"""
Phase 2 — manual centrality calculation (bypass igraph eigenvector bug)
"""
import sys, time
import numpy as np
import pandas as pd
import networkx as nx
import celloracle as co

t0 = time.time()
def stamp(m): print(f"[{time.time()-t0:6.1f}s] {m}", flush=True)

BASE = "."

stamp(f"celloracle {co.__version__}")
oracle = co.load_hdf5(f"{BASE}/r1r9.celloracle.oracle")
links  = co.load_hdf5(f"{BASE}/r1r9.celloracle.links")
stamp("Loaded oracle + links")

# Re-filter (overwritten from previous run probably, but redo)
links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)


# ---- Manual centrality via NetworkX --------------------------
def compute_centrality(edge_df):
    """edge_df: source, target, coef_mean, coef_abs, p"""
    G = nx.DiGraph()
    for _, row in edge_df.iterrows():
        G.add_edge(row["source"], row["target"],
                   weight=abs(row["coef_mean"]))
    nodes = list(G.nodes())

    deg_in  = dict(G.in_degree())
    deg_out = dict(G.out_degree())
    deg_all = {n: deg_in.get(n,0) + deg_out.get(n,0) for n in nodes}

    # Eigenvector on undirected weighted
    Gu = G.to_undirected()
    try:
        eig = nx.eigenvector_centrality_numpy(Gu, weight="weight")
    except Exception:
        eig = {n: np.nan for n in nodes}

    # Betweenness (unweighted for speed)
    bet = nx.betweenness_centrality(G, k=None, normalized=True)

    return pd.DataFrame({
        "gene": nodes,
        "degree_centrality_in":  [deg_in.get(n,0) for n in nodes],
        "degree_centrality_out": [deg_out.get(n,0) for n in nodes],
        "degree_centrality_all": [deg_all[n] for n in nodes],
        "eigenvector_centrality":[eig.get(n,np.nan) for n in nodes],
        "betweenness_centrality":[bet.get(n,np.nan) for n in nodes],
    })

stamp("Computing centrality per cluster (NetworkX) ...")
scores = []
for cl, df in links.filtered_links.items():
    s = compute_centrality(df)
    s["cluster"] = cl
    scores.append(s)
    stamp(f"  {cl}: {len(s)} nodes")
all_scores = pd.concat(scores, ignore_index=True)
all_scores.to_csv(f"{BASE}/phase2_centrality.csv", index=False)
stamp("Saved centrality table")


# ---- Top hubs per cluster ----------------------------------
TOP_N = 20
print("\n========= Top hubs per cluster (eigenvector) =========\n")
for cl in sorted(links.filtered_links.keys()):
    sub = all_scores[all_scores["cluster"]==cl].copy()
    top = sub.sort_values("eigenvector_centrality", ascending=False).head(TOP_N)
    print(f"\n--- {cl} (top-{TOP_N} by eigenvector) ---")
    print(top[["gene","eigenvector_centrality",
               "degree_centrality_all","betweenness_centrality"]]
          .to_string(index=False))

print("\n========= Top hubs per cluster (degree_all) =========\n")
for cl in sorted(links.filtered_links.keys()):
    sub = all_scores[all_scores["cluster"]==cl].copy()
    top = sub.sort_values("degree_centrality_all", ascending=False).head(TOP_N)
    print(f"\n--- {cl} (top-{TOP_N} by degree_all) ---")
    print(top[["gene","degree_centrality_all",
               "degree_centrality_out","eigenvector_centrality"]]
          .to_string(index=False))


# ---- Lipid TF cross-check ----------------------------------
lipid_tfs = ["Srebf1","Srebf2","Pparg","Ppara","Ppard",
             "Mlxipl","Mlx","Nr1h2","Nr1h3","Nr2f2",
             "Insig1","Insig2","Cebpa","Cebpb","Klf2","Klf4",
             "Egr1","Stat1","Stat3","Foxo1","Foxo3",
             "Atf4","Xbp1","Ncoa3","Ppargc1a"]

lipid_sub = all_scores[all_scores["gene"].isin(lipid_tfs)]
piv_eig = lipid_sub.pivot(index="gene", columns="cluster",
                          values="eigenvector_centrality").round(3)
piv_deg = lipid_sub.pivot(index="gene", columns="cluster",
                          values="degree_centrality_all").astype("Int64")

print("\n========= Lipid TFs — eigenvector centrality =========")
print(piv_eig)
print("\n========= Lipid TFs — degree_all (out-degree as TF) ===")
piv_out = lipid_sub.pivot(index="gene", columns="cluster",
                          values="degree_centrality_out").astype("Int64")
print(piv_out)

piv_eig.to_csv(f"{BASE}/phase2_lipid_eigenvector.csv")
piv_out.to_csv(f"{BASE}/phase2_lipid_degree_out.csv")

stamp(f"DONE — total {time.time()-t0:.1f}s")
