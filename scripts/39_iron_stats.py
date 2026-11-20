#!/usr/bin/env python
"""
39_iron_stats.py
Iron dose response (Results, "Complementary antagonistic mechanisms"; Fig 3B).
Each yeast-pathogen pairing is one observation: did mean %MI fall between 5 and
20 mg/l FeCl3? Two-sided sign test (exact binomial on the direction of change) and
Wilcoxon signed-rank test, with B. cinerea and the two Fusarium pathogens analysed
separately.
Input:  analysis/data/phenotypes_long.csv (assay == "iron")
Output: analysis/data/iron_stats.csv (+ printed summary)
"""
import numpy as np
import pandas as pd
from scipy import stats

DATA = "analysis/data"
d = pd.read_csv(f"{DATA}/phenotypes_long.csv")
d = d[d.assay == "iron"]
w = d.pivot_table(index=["isolate", "pathogen"], columns="condition",
                  values="pct_MI").reset_index()
w["change"] = w["20 mg/L Fe"] - w["5 mg/L Fe"]
w["group"] = np.where(w.pathogen == "B. cinerea", "B. cinerea", "Fusarium")

rows = []
for g, x in w.groupby("group"):
    ch = x.change.values
    nz = ch[ch != 0]
    n_down = int((nz < 0).sum())
    sign_p = stats.binomtest(n_down, len(nz), 0.5, alternative="two-sided").pvalue
    wil_p = stats.wilcoxon(ch, alternative="two-sided").pvalue
    rows.append(dict(group=g, n=len(ch), n_decreased=n_down, n_increased=int((nz > 0).sum()),
                     median_change_pp=round(float(np.median(ch)), 1),
                     min_change_pp=round(float(ch.min()), 1),
                     max_change_pp=round(float(ch.max()), 1),
                     sign_test_p=round(sign_p, 3), wilcoxon_p=round(wil_p, 3)))
    print(f"\n{g}:")
    print(x[["isolate", "pathogen", "5 mg/L Fe", "20 mg/L Fe", "change"]]
          .round(1).to_string(index=False))

out = pd.DataFrame(rows)
print("\n" + out.to_string(index=False))
out.to_csv(f"{DATA}/iron_stats.csv", index=False)
