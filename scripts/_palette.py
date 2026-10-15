"""Shared species colour palette + label helpers for all figures.
One colour per species across every figure (dataviz: colour follows the entity).
Okabe-Ito colourblind-safe hues in a fixed order."""

# fixed species order (by overall prevalence / taxonomic grouping)
SP_ORDER = ["Hanseniaspora uvarum","Sporobolomyces roseus","Rhodotorula babjevae",
            "Pichia sp.","Metschnikowia sp.","Vishniacozyma carnescens",
            "Sporobolomyces sp.","Naganishia vaughanmartiniae","Filobasidium stepposum",
            "Hannaella sp."]
_OI = ["#0072B2","#E69F00","#009E73","#CC79A7","#D55E00","#56B4E9","#F0E442",
       "#999999","#8c564b","#17becf"]
SP_COL = {s:_OI[i] for i,s in enumerate(SP_ORDER)}
SP_COL["Vishniacozyma victoriae"] = "#17becf"   # sister of carnescens (cyan vs sky-blue)
SP_COL["Filobasidium sp."] = "#8c564b"
SP_COL["Fusarium venenatum (outgroup)"] = "#000000"
SP_COL["Other"] = "#bbbbbb"

def spc(s):
    """Species colour; also accepts the cf. display names (e.g. 'Rhodotorula aff. babjevae')."""
    if s in SP_COL:
        return SP_COL[s]
    for k, v in CF_NAME.items():
        if v == s and k in SP_COL:
            return SP_COL[k]
    return "#bbbbbb"

# Display names applied consistently across every figure. The unresolved Rhodotorula
# lineage uses open nomenclature (aff.); Metschnikowia is assigned to M. pulcherrima.
CF_NAME = {
    "Rhodotorula babjevae": "Rhodotorula aff. babjevae",
    "Metschnikowia sp.": "Metschnikowia pulcherrima",
    "Metschnikowia sp": "Metschnikowia pulcherrima",
    "Metschnikowia ziziphicola": "Metschnikowia pulcherrima",
    "Pichia sp.": "Pichia kluyveri",
    "Pichia sp": "Pichia kluyveri",
    "Sporobolomyces sp.": "Sporobolomyces sp.",
}


def disp(s):
    """Canonical display name for a species label."""
    return CF_NAME.get(str(s).strip(), s)


def it(s):
    """Italicise a binomial in mathtext, leaving sp., cf., aff. and outgroup upright."""
    s = disp(s)
    if "outgroup" in s:
        return s
    parts = []
    for w in s.split(" "):
        if w in ("sp.", "sp", "cf.", "aff."):
            parts.append(w)
        else:
            parts.append(f"$\\it{{{w}}}$")
    return " ".join(parts)
