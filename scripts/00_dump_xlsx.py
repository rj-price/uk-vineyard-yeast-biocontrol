import openpyxl, glob, os, csv
files = sorted(glob.glob("gideon_mphil/*.xlsx") + glob.glob("genomics/*.xlsx"))
outdir = "analysis/data/raw_csv"
os.makedirs(outdir, exist_ok=True)
for f in files:
    base = os.path.splitext(os.path.basename(f))[0].replace(" ","_").replace("&","and")
    wb = openpyxl.load_workbook(f, data_only=True)
    for ws in wb.worksheets:
        safe = ws.title.replace(" ","_").replace(".","").replace("/","-")
        out = os.path.join(outdir, f"{base}__{safe}.csv")
        with open(out,"w",newline="") as fh:
            w = csv.writer(fh)
            for row in ws.iter_rows(values_only=True):
                w.writerow(["" if v is None else v for v in row])
        print("wrote", out)
