#!/usr/bin/env python
"""
37_score_sensitivity.py
Robustness of the Fig 6 composite prioritisation score (reviewer point 9).

(1) Weights: every weight vector with each of the five traits in {0, 1, 2}
    (3^5 - 1 = 242 schemes, all-zero excluded).
(2) Replicates: 2,000 bootstrap draws resampling the three replicate plates of each
    isolate x assay x pathogen x condition, at the published weights.
The score definition matches 02_stats_figures.py. Enzyme breadth has no replicate
plates in the tidy table and is held fixed in (2).

-> analysis/data/score_sensitivity.csv  (per isolate: published rank, % top-1/top-3
   across weight schemes and bootstrap draws, bootstrap 95% rank interval)
"""
import itertools
import numpy as np, pandas as pd

PROJ = "/mnt/shared/scratch/jnprice/private/biocontrol_yeast/analysis"
rng = np.random.default_rng(7)
NBOOT = 2000
TRAITS = ["broad_spectrum_MI", "MI_Bcinerea", "VOC_mean", "Fe_comp_mean", "enzyme_breadth"]
W_PUB = np.array([2, 1, 1, 1, 1])

long = pd.read_csv(f"{PROJ}/data/phenotypes_long.csv")
wide = pd.read_csv(f"{PROJ}/data/phenotype_matrix.csv").set_index("isolate")
iso = wide.index[wide.MI_Bcinerea.notna()]
enz = (wide.loc[iso, [c for c in wide if c.startswith("enz_")]] > 1.05).sum(axis=1)

# control diameter per row, recovered from the published mean and %MI
long = long[long.isolate.isin(iso)].copy()
long["Dc"] = np.where(long.pct_MI < 100, long.mean_growth / (1 - long.pct_MI / 100), np.nan)
long["Dc"] = long.Dc.fillna(85.0)
R = long[["R1", "R2", "R3"]].values


def traits(growth):
    """growth: per-row mean pathogen diameter -> trait table indexed by isolate."""
    d = long.assign(mi=(long.Dc - growth) / long.Dc * 100)
    piv = lambda a: d[d.assay == a].pivot_table(index="isolate", columns=["pathogen", "condition"],
                                                 values="mi", dropna=False)
    dual = d[d.assay == "dual"].pivot_table(index="isolate", columns="pathogen", values="mi")
    t = pd.DataFrame(index=iso)
    t["MI_Bcinerea"] = dual["B. cinerea"]
    t["broad_spectrum_MI"] = dual[["B. cinerea", "F. oxysporum", "F. graminearum"]].mean(axis=1)
    t["VOC_mean"] = d[d.assay == "VOC"].groupby("isolate").mi.mean()
    t["Fe_comp_mean"] = d[d.assay == "iron"].groupby("isolate").mi.mean()
    t["enzyme_breadth"] = enz
    return t


def score(t, w):
    z = (t - t.mean()) / t.std()
    return (z.fillna(0)[TRAITS] * w).sum(axis=1)


def ranks(s):
    return s.rank(ascending=False, method="min")


base = traits(long.mean_growth.values)
pub = ranks(score(base, W_PUB))

# (1) weight schemes
W = [np.array(w) for w in itertools.product([0, 1, 2], repeat=5) if any(w)]
wr = pd.DataFrame({i: ranks(score(base, w)) for i, w in enumerate(W)})

# (2) replicate bootstrap
br = {}
for b in range(NBOOT):
    pick = rng.integers(0, 3, size=R.shape)
    g = np.take_along_axis(R, pick, axis=1).mean(axis=1)
    br[b] = ranks(score(traits(g), W_PUB))
br = pd.DataFrame(br)

out = pd.DataFrame({
    "species": wide.loc[iso, "species"],
    "published_rank": pub.astype(int),
    "weights_pct_top1": (wr == 1).mean(axis=1).mul(100).round(1),
    "weights_pct_top3": (wr <= 3).mean(axis=1).mul(100).round(1),
    "weights_median_rank": wr.median(axis=1),
    "boot_pct_top1": (br == 1).mean(axis=1).mul(100).round(1),
    "boot_pct_top3": (br <= 3).mean(axis=1).mul(100).round(1),
    "boot_rank_2.5": br.quantile(0.025, axis=1),
    "boot_rank_97.5": br.quantile(0.975, axis=1),
}).sort_values("published_rank")
out.to_csv(f"{PROJ}/data/score_sensitivity.csv")
print(f"{len(W)} weight schemes, {NBOOT} bootstrap draws")
print(out.head(10).to_string())
