#!/usr/bin/env python
"""
02_stats_figures.py
Statistics + publication figures for the yeast biocontrol screen.
All figures -> analysis/figures/ ; stats tables -> analysis/data/.
"""
import os, itertools, warnings
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
from statsmodels.stats.multicomp import pairwise_tukeyhsd
warnings.filterwarnings("ignore")

DATA="analysis/data"; FIG="analysis/figures"; os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,
                     "axes.grid":True,"grid.alpha":0.25,"grid.linewidth":0.6,
                     "figure.dpi":140,"savefig.dpi":300,"savefig.bbox":"tight",
                     "font.family":"DejaVu Sans"})

import seaborn as sns
import sys; sys.path.insert(0, os.path.dirname(__file__))
from _palette import SP_COL as sp_col, spc, it

master = pd.read_csv(f"{DATA}/isolate_master.csv")
long   = pd.read_csv(f"{DATA}/phenotypes_long.csv")
wide   = pd.read_csv(f"{DATA}/phenotype_matrix.csv")

# species present among phenotyped isolates, in the shared fixed order
from _palette import SP_ORDER
sp_order = [s for s in SP_ORDER if s in set(master["species"].dropna())
            and s in {"Hanseniaspora uvarum","Rhodotorula babjevae","Pichia sp.",
                      "Metschnikowia sp.","Vishniacozyma carnescens"}]

# ============================================================================
# FIGURE 1 — Antagonism (dual-plate %MI) by pathogen, points coloured by species
# ============================================================================
dual = long[long.assay=="dual"].copy()
dual = dual.merge(master[["isolate","species"]], on="isolate", how="left")
paths = ["B. cinerea","F. oxysporum","F. graminearum"]
fig, ax = plt.subplots(figsize=(7.2,4.6))
rng = np.random.default_rng(0)
for i,p in enumerate(paths):
    sub = dual[dual.pathogen==p]
    ax.boxplot(sub.pct_MI, positions=[i], widths=0.55, showfliers=False,
               patch_artist=True, boxprops=dict(facecolor="#eeeeee",edgecolor="#555"),
               medianprops=dict(color="#222",lw=1.6), whiskerprops=dict(color="#777"),
               capprops=dict(color="#777"))
    x = i + (rng.random(len(sub))-0.5)*0.28
    ax.scatter(x, sub.pct_MI, s=34, c=[spc(s) for s in sub.species],
               edgecolor="white", linewidth=0.6, zorder=3)
ax.set_xticks(range(len(paths))); ax.set_xticklabels([it(p) for p in paths])
ax.set_ylabel("Mycelial growth inhibition (%)")
ax.set_title("Dual-culture antagonism varies by pathogen and yeast species", fontsize=11, loc="left")
ax.axhline(0, color="#bbb", lw=0.8)
handles=[Patch(facecolor=sp_col[s],label=it(s)) for s in sp_order]
ax.legend(handles=handles, title="Species", frameon=False, fontsize=8, title_fontsize=8,
          loc="upper right")
plt.savefig(f"{FIG}/fig1_antagonism_by_pathogen.png"); plt.close()

# ============================================================================
# FIGURE 2 — clustered heatmap of all phenotypic traits (isolate x trait)
# ============================================================================
trait_cols = [c for c in wide.columns if c.startswith(("MI_","VOC_","Fe_","enz_"))]
H = wide.set_index("isolate")[trait_cols].apply(pd.to_numeric, errors="coerce")
H = H.loc[H.notna().sum(axis=1) >= 3]                 # keep well-characterised isolates
# z-score each trait for cross-trait comparability
Hz = (H - H.mean()) / H.replace(0,np.nan).std().replace(0,1)
Hz = Hz.fillna(0)
row_species = wide.set_index("isolate").loc[Hz.index,"species"].fillna("Other")
rc = [spc(s) for s in row_species]
pretty = {c: c.replace("MI_","Dual: ").replace("VOC_","VOC: ").replace("Fe_","Fe: ")
             .replace("enz_","Enz: ").replace("Bcinerea","B. cinerea")
             .replace("Fgraminearum","F. graminearum").replace("Foxysporum","F. oxysporum")
             .replace("_5mgLFe"," (5)").replace("_20mgLFe"," (20)") for c in Hz.columns}
g = sns.clustermap(Hz.rename(columns=pretty), cmap="RdBu_r", center=0, vmin=-2.2, vmax=2.2,
                   row_colors=rc, figsize=(11,9.5), linewidths=.4, linecolor="white",
                   cbar_kws={"label":"z-score (per trait)"}, dendrogram_ratio=(.12,.10),
                   col_cluster=True, row_cluster=True)
g.ax_heatmap.set_xlabel(""); g.ax_heatmap.set_ylabel("Isolate")
for lab in g.ax_heatmap.get_xticklabels(): lab.set_rotation(60); lab.set_ha("right"); lab.set_fontsize(8)
handles=[Patch(facecolor=sp_col[s],label=it(s)) for s in sp_order]
g.ax_heatmap.legend(handles=handles, title="Species", bbox_to_anchor=(1.28,1.0),
                    loc="upper left", frameon=False, fontsize=8, title_fontsize=8)
g.fig.suptitle("Phenotypic trait profiles cluster yeast isolates by functional type",
               x=0.5,y=1.02,fontsize=12)
g.savefig(f"{FIG}/fig2_trait_heatmap.png"); plt.close()

# ============================================================================
# FIGURE 3 — PCA of isolates in phenotype space
# ============================================================================
X = H.copy()
X = X.loc[:, X.notna().sum() >= X.shape[0]*0.6]       # traits present for most isolates
X = X.dropna(axis=0, how="any") if X.isna().sum().sum() else X
X = X.fillna(X.mean())
Xs = StandardScaler().fit_transform(X)
pca = PCA(n_components=3).fit(Xs); PC = pca.transform(Xs)
ev = pca.explained_variance_ratio_*100
pcdf = pd.DataFrame(PC[:,:2], columns=["PC1","PC2"], index=X.index)
pcdf["species"] = wide.set_index("isolate").loc[pcdf.index,"species"].fillna("Other")
fig, ax = plt.subplots(figsize=(6.6,5.4))
for s in sp_order+["Other"]:
    sub = pcdf[pcdf.species==s]
    if len(sub): ax.scatter(sub.PC1, sub.PC2, s=70, color=spc(s), edgecolor="white",
                            linewidth=0.8, label=it(s), zorder=3)
for iso,r in pcdf.iterrows():
    ax.annotate(iso,(r.PC1,r.PC2), fontsize=6.5, color="#444",
                xytext=(3,3), textcoords="offset points")
# loadings
ld = pca.components_[:2].T * np.sqrt(pca.explained_variance_[:2])
sc = 3.2/np.abs(ld).max()
for j,c in enumerate(X.columns):
    ax.arrow(0,0,ld[j,0]*sc,ld[j,1]*sc, color="#999", alpha=0.5, head_width=0.06, lw=0.6, zorder=1)
ax.set_xlabel(f"PC1 ({ev[0]:.0f}%)"); ax.set_ylabel(f"PC2 ({ev[1]:.0f}%)")
ax.set_title("PCA separates strong broad-spectrum antagonists", fontsize=11, loc="left")
ax.legend(frameon=False, fontsize=8, title="Species", title_fontsize=8, loc="best")
ax.axhline(0,color="#ddd",lw=0.7); ax.axvline(0,color="#ddd",lw=0.7)
plt.savefig(f"{FIG}/fig3_pca.png"); plt.close()

# ============================================================================
# FIGURE 4 — trait-trait correlation (diverging), Spearman
# ============================================================================
corr = H.corr(method="spearman")
fig, ax = plt.subplots(figsize=(9,7.5))
im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1)
ax.set_xticks(range(len(corr))); ax.set_yticks(range(len(corr)))
ax.set_xticklabels([pretty[c] for c in corr.columns], rotation=60, ha="right", fontsize=7)
ax.set_yticklabels([pretty[c] for c in corr.index], fontsize=7)
for i in range(len(corr)):
    for j in range(len(corr)):
        v=corr.iloc[i,j]
        if not np.isnan(v): ax.text(j,i,f"{v:.1f}",ha="center",va="center",
                                    fontsize=6,color="white" if abs(v)>0.55 else "#333")
cb=plt.colorbar(im,fraction=0.046,pad=0.04); cb.set_label("Spearman ρ")
ax.set_title("Correlation among biocontrol traits", fontsize=11, loc="left")
plt.savefig(f"{FIG}/fig4_trait_correlation.png"); plt.close()

# ============================================================================
# STATS 1 — one-way ANOVA + Tukey on dual-plate B. cinerea %MI across isolates
#          (replicate-level; %MI computed per rep vs control mean 85)
# ============================================================================
def rep_MI(df):  # long rows already have R1..R3 growth; MI per rep
    out=[]
    for _,r in df.iterrows():
        for rp in ["R1","R2","R3"]:
            out.append((r.isolate, (85-r[rp])/85*100))
    return pd.DataFrame(out, columns=["isolate","MI"])
statlines=[]
for p in paths:
    d = rep_MI(long[(long.assay=="dual")&(long.pathogen==p)])
    groups=[g.MI.values for _,g in d.groupby("isolate")]
    F,pv = stats.f_oneway(*groups)
    statlines.append(f"Dual {p}: one-way ANOVA across isolates F={F:.1f}, p={pv:.2e}")
# Tukey for B. cinerea (most active)
dbc = rep_MI(long[(long.assay=="dual")&(long.pathogen=="B. cinerea")])
tuk = pairwise_tukeyhsd(dbc.MI, dbc.isolate, alpha=0.05)
tuk_df = pd.DataFrame(tuk.summary().data[1:], columns=tuk.summary().data[0])
tuk_df.to_csv(f"{DATA}/tukey_dual_Bcinerea.csv", index=False)

# STATS 2 — %MI by species (dual, B. cinerea), Kruskal-Wallis
dsp = dual[dual.pathogen=="B. cinerea"].dropna(subset=["species"])
kw = stats.kruskal(*[g.pct_MI.values for _,g in dsp.groupby("species") if len(g)>=2])
statlines.append(f"Dual B. cinerea %MI by species: Kruskal-Wallis H={kw.statistic:.2f}, p={kw.pvalue:.3f}")

# ============================================================================
# CANDIDATE RANKING FRAMEWORK
# ============================================================================
rank = wide.set_index("isolate").copy()
rank["broad_spectrum_MI"] = rank[["MI_Bcinerea","MI_Foxysporum","MI_Fgraminearum"]].mean(axis=1)
rank["VOC_mean"] = rank[[c for c in rank if c.startswith("VOC_")]].mean(axis=1)
rank["Fe_comp_mean"] = rank[[c for c in rank if c.startswith("Fe_")]].mean(axis=1)
enzc=[c for c in rank if c.startswith("enz_")]
rank["enzyme_breadth"] = (rank[enzc] > 1.05).sum(axis=1)          # # enzymes active
def z(s):
    s=pd.to_numeric(s,errors="coerce"); return (s-s.mean())/s.std()
rank["score"] = (z(rank.broad_spectrum_MI).fillna(0)*2 + z(rank.MI_Bcinerea).fillna(0)
                 + z(rank.VOC_mean).fillna(0) + z(rank.Fe_comp_mean).fillna(0)
                 + z(rank.enzyme_breadth).fillna(0))
rank = rank.join(master.set_index("isolate")[["species","seq_id","genome_code"]], rsuffix="_m")
cand = rank[["species","seq_id","MI_Bcinerea","MI_Foxysporum","MI_Fgraminearum",
             "broad_spectrum_MI","VOC_mean","Fe_comp_mean","enzyme_breadth","score"]]
cand = cand.dropna(subset=["MI_Bcinerea"]).sort_values("score", ascending=False).round(2)
cand.to_csv(f"{DATA}/candidate_ranking.csv")

# ranking figure (top candidates)
top = cand.head(12).iloc[::-1]
fig, ax = plt.subplots(figsize=(7.2,5))
colors=[spc(s) for s in top.species]
ax.barh(range(len(top)), top.score, color=colors, edgecolor="white")
ax.set_yticks(range(len(top))); ax.set_yticklabels(top.index)
ax.set_xlabel("Composite biocontrol score (z-weighted)")
ax.set_title("Candidate ranking: top broad-spectrum yeast isolates", fontsize=11, loc="left")
handles=[Patch(facecolor=sp_col[s],label=it(s)) for s in sp_order if s in set(top.species)]
ax.legend(handles=handles, frameon=False, fontsize=8, title="Species", title_fontsize=8, loc="lower right")
plt.savefig(f"{FIG}/fig5_candidate_ranking.png"); plt.close()

with open(f"{DATA}/stats_summary.txt","w") as fh:
    fh.write("\n".join(statlines)+"\n")
print("\n".join(statlines))
print("\nTop candidates:\n", cand.head(8).to_string())
print("\nFigures written to", FIG)
