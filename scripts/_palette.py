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
    return SP_COL.get(s, "#bbbbbb")

def it(s):
    """Italicise a binomial in mathtext (keeping spaces); leave 'sp.'/outgroup upright."""
    if "sp." in s or "outgroup" in s:
        return s
    return f"$\\it{{{s.replace(' ', chr(92)+' ')}}}$"
