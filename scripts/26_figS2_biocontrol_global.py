#!/usr/bin/env python
"""
26_figS2_biocontrol_global.py -> figures/manuscript/FigS2_biocontrol_global.{svg,png}
Supplementary figure — biocontrol-relevant gene content is indistinguishable between the
Kent isolates and their global conspecifics, for both lead lineages. This is the
gene-content complement to the ANI/phylogeny (Fig S1): the lineages are cosmopolitan AND
carry the same biocontrol repertoire wherever they come from — consistent with the paper's
central finding that Botrytis antagonism is not readable from gene content.

Parses, per mined genome (analysis/mining/results/<iso>/), exactly as 08/10:
  dbcan/overview.tsv       -> CAZyme families (>=2 tools) -> enzyme classes (total & secreted)
  <iso>.annotations.txt    -> MEROPS proteases (total & secreted), CAZyme per-transcript
  <iso>.secreted_ids.txt   -> SignalP Sec/SPI transcript IDs
  antismash/<iso>.json     -> secondary-metabolite cluster count
"""
import os, re, sys, json, glob, collections
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
sys.path.insert(0, os.path.dirname(__file__))
import _figstyle as fs

MINE = "analysis/mining"
OUT = "analysis/figures/manuscript"
DATA = "analysis/data"

# genus -> (Kent isolates, [(global label, results-id, origin)])
GROUPS = {
    "Rhodotorula sp. (aff. babjevae)": (
        ["UY1", "TY27", "TY28", "TY29"],
        [("JJ10.1", "RHOjj101", "China"), ("DBVPG 8058", "RHOdbvpg", "coll.")]),
    "Hanseniaspora uvarum": (
        ["RY14", "RY3", "RY4", "RY7", "RY25"],
        [("UKR", "HUVukr", "Ukraine"), ("AUS", "HUVaus", "Australia"),
         ("CHN", "HUVchn", "China")]),
}

CAZY_CLASS = {
    "chitinase\n(GH18/19)": ["GH18", "GH19"],
    "β-1,3-glucanase": ["GH16", "GH17", "GH55", "GH64", "GH81", "GH128", "GH152"],
    "β-1,6-glucanase": ["GH30"],
    "cellulase": ["GH5", "GH6", "GH7", "GH12", "GH45"],
    "amylase": ["GH13", "GH15"],
    "pectinase": ["GH28", "PL1", "PL3", "PL4"],
}
FAM2CLASS = {f: c for c, fl in CAZY_CLASS.items() for f in fl}


def fams_from(cell):
    return re.findall(r"(GH\d+|GT\d+|CBM\d+|CE\d+|PL\d+|AA\d+)", str(cell))


def inventory(iso):
    d = f"{MINE}/results/{iso}"
    ann = pd.read_csv(f"{d}/{iso}.annotations.txt", sep="\t", dtype=str).fillna("")
    sec = set(l.strip() for l in open(f"{d}/{iso}.secreted_ids.txt") if l.strip())
    idcol = "TranscriptID" if "TranscriptID" in ann.columns else ann.columns[1]
    ann["sec"] = ann[idcol].isin(sec)
    rec = {"secretome": int(ann["sec"].sum()), "n_genes": len(ann)}
    # CAZyme classes (total + secreted)
    for cls in CAZY_CLASS:
        rec[cls] = 0
        rec[cls + "|sec"] = 0
    # same CAZyme call as Table S6 (08_mining_summary.py): dbCAN genes supported by
    # >= 2 tools, family taken from dbCAN_hmm, then DIAMOND, then dbCAN_sub
    ov = pd.read_csv(f"{d}/dbcan/overview.tsv", sep="\t")
    ov.columns = [c.strip() for c in ov.columns]
    ov = ov[pd.to_numeric(ov["#ofTools"], errors="coerce").fillna(0) >= 2]
    for _, r in ov.iterrows():
        fam = next((f[0] for f in (fams_from(r[c]) for c in ("dbCAN_hmm", "DIAMOND", "dbCAN_sub")) if f), None)
        c = FAM2CLASS.get(fam)
        if c:
            rec[c] += 1
            if r["Gene ID"] in sec:
                rec[c + "|sec"] += 1
    # proteases (MEROPS col)
    pcol = "Protease" if "Protease" in ann.columns else None
    if pcol:
        isp = ann[pcol] != ""
        rec["protease"] = int(isp.sum())
        rec["protease|sec"] = int((isp & ann["sec"]).sum())
    # antiSMASH clusters
    js = glob.glob(f"{d}/antismash/*.json")
    nsm = 0
    if js:
        data = json.load(open(js[0]))
        for r in data.get("records", []):
            nsm += len(r.get("areas", []))
    rec["SM clusters"] = nsm
    return rec


def build():
    rows = []
    for genus, (kent, glob_) in GROUPS.items():
        for iso in kent:
            rows.append(dict(genus=genus, label=iso, origin="Kent, UK",
                             grp="Kent", **inventory(iso)))
        for lab, gid, origin in glob_:
            rows.append(dict(genus=genus, label=f"{lab}", origin=origin,
                             grp="global", **inventory(gid)))
    return pd.DataFrame(rows)


DISPLAY = ["secretome", "protease", "protease|sec", "SM clusters",
           "chitinase\n(GH18/19)", "β-1,3-glucanase", "β-1,3-glucanase|sec",
           "β-1,6-glucanase", "amylase", "pectinase"]


def panel(ax, df, genus, cax=None):
    sub = df[df.genus == genus].reset_index(drop=True)
    cols = [c for c in DISPLAY if c in sub.columns]
    M = sub[cols].astype(float).values
    # colour = deviation from the column median on a fixed 0-30 % scale (as Fig 4D), so
    # a one- or two-gene difference in a large count stays pale; raw counts annotated
    Z = np.zeros_like(M)
    for j in range(M.shape[1]):
        med = np.nanmedian(M[:, j])
        dev = np.abs(M[:, j] - med)
        # a one-gene difference is the annotation noise floor: left unshaded
        Z[:, j] = np.where(dev <= 1, 0.0, dev / med if med else 0.0)
    VMAX = 0.30
    im = ax.imshow(Z, aspect="auto", cmap="Blues", vmin=0, vmax=VMAX)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c.replace("|sec", " (secr.)").replace("\n", " ") for c in cols],
                       fontsize=5.8, rotation=35, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels([f"{r.label}  ({r.origin})" for r in sub.itertuples()], fontsize=6)
    for i in range(len(sub)):
        for j in range(len(cols)):
            ax.text(j, i, f"{int(M[i,j])}", ha="center", va="center", fontsize=5.4,
                    color="white" if Z[i, j] > 0.6 * VMAX else fs.INK)
    nk = (sub.grp == "Kent").sum()
    ax.axhline(nk - 0.5, color="#CC3311", lw=1.8)
    yt = ax.get_yaxis_transform()   # x in axes fraction: sit left of the isolate names
    ax.text(-0.25, (nk - 1) / 2, "Kent", rotation=90, va="center", ha="right",
            fontsize=6, color="#CC3311", fontweight="bold", transform=yt)
    ax.text(-0.25, nk + (len(sub) - nk - 1) / 2, "global", rotation=90, va="center",
            ha="right", fontsize=6, color=fs.MUTED, fontweight="bold", transform=yt)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(fs.patho_it(genus) if genus.startswith("Hansen") else
                 "$\\it{Rhodotorula}$ sp. (aff. $\\it{babjevae}$)", fontsize=8, pad=16)
    if cax is not None:
        cb = plt.colorbar(im, cax=cax, ticks=[0, 0.15, 0.30])
        cb.ax.set_yticklabels(["0", "15%", "30%"]); cb.ax.tick_params(labelsize=5.5, length=2)
        cb.set_label("deviation from column median", fontsize=5.5)


def stats_note(df, genus):
    # "secretome" is total SignalP count (tracks total gene number, assembly-sensitive);
    # kept separate from the specific biocontrol enzyme/protease/SM classes.
    sub = df[df.genus == genus]
    k = sub[sub.grp == "Kent"]; g = sub[sub.grp == "global"]
    bits = []
    for c in DISPLAY:
        if c not in sub.columns or c == "secretome":
            continue
        try:
            p = mannwhitneyu(k[c], g[c], alternative="two-sided").pvalue
        except ValueError:
            p = np.nan
        if p < 0.05:
            bits.append(f"{c.splitlines()[0]} p={p:.3f}")
    # State the smallest P value these group sizes can reach, so a null is not read as
    # evidence of absence. With n_Kent vs n_global genomes, a two-sided Mann-Whitney U
    # bottoms out at the value below even under perfect separation.
    pfloor = mannwhitneyu(list(range(len(k))),
                          list(range(100, 100 + len(g))),
                          alternative="two-sided").pvalue
    note = ("no biocontrol category differs (Mann–Whitney; n=%d vs %d, "
            "smallest attainable p=%.3f)" % (len(k), len(g), pfloor)
            if not bits else "differs: " + "; ".join(bits) +
            " (n=%d vs %d, smallest attainable p=%.3f)" % (len(k), len(g), pfloor))
    # flag the secretome-size proxy separately if it separates the groups
    try:
        ps = mannwhitneyu(k["secretome"], g["secretome"], alternative="two-sided").pvalue
    except ValueError:
        ps = np.nan
    if ps < 0.06:
        note += ("; total secretome size p=%.3f (tracks total gene count, "
                 "not a specific mechanism)" % ps)
    return note


def main():
    os.makedirs(OUT, exist_ok=True)
    fs.apply()
    df = build()
    df.to_csv(f"{DATA}/biocontrol_global_inventory.csv", index=False)

    fig = plt.figure(figsize=(7.4, 7.4))
    axR = fig.add_axes([0.22, 0.62, 0.60, 0.25])
    axH = fig.add_axes([0.22, 0.16, 0.60, 0.27])
    caxR = fig.add_axes([0.85, 0.64, 0.015, 0.20])
    panel(axR, df, "Rhodotorula sp. (aff. babjevae)", caxR)
    panel(axH, df, "Hanseniaspora uvarum")

    fig.text(0.02, 0.93, "A", fontsize=12, fontweight="bold", va="top")
    fig.text(0.02, 0.49, "B", fontsize=12, fontweight="bold", va="top")
    for genus, y in [("Rhodotorula sp. (aff. babjevae)", 0.515),
                     ("Hanseniaspora uvarum", 0.055)]:
        fig.text(0.22, y, stats_note(df, genus), fontsize=6, color=fs.MUTED, va="top")
    fig.suptitle("Biocontrol gene content: Kent isolates vs global conspecifics",
                 fontsize=9, y=0.955, x=0.5)
    fs.save(fig, f"{OUT}/FigS2_biocontrol_global")
    print("wrote FigS2_biocontrol_global")
    print(df[["genus", "label", "grp"] + [c for c in DISPLAY if c in df.columns]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
