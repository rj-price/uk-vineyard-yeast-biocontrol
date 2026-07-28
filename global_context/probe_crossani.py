import pandas as pd, os
GC = "analysis/global_context"
names = pd.read_csv(f"{GC}/rhodo_names.csv")
p2n = dict(zip(names.path.map(os.path.abspath), names.name))
a = pd.read_csv(f"{GC}/rhodo_ani.tsv", sep="\t", header=None,
                names=["q", "r", "ani", "frag", "tot"])
a["qn"] = a.q.map(lambda p: p2n.get(os.path.abspath(p)))
a["rn"] = a.r.map(lambda p: p2n.get(os.path.abspath(p)))

lab = {"GCA_920104245.2": "babjevae-TYPE (CBS7808)",
       "GCA_050403105.1": "babjevae JCM9279",
       "GCA_920103735.2": "'babjevae' DBVPG8058",
       "GCA_026119225.1": "'graminis' JJ10.1 China",
       "GCA_937859235.1": "'graminis' MAG Australia",
       "GCA_001329695.1": "graminis-REF WP1",
       "GCA_050403165.1": "graminis JCM3775",
       "UY1": "Kent UY1"}
keys = list(lab)


def ani(x, y):
    s = a[(a.qn == x) & (a.rn == y)]
    return f"{s.ani.iloc[0]:.1f}" if len(s) else "  <80"


print(f"{'':26}" + "".join(f"{lab[k][:11]:>12}" for k in keys))
for x in keys:
    print(f"{lab[x]:26}" + "".join(f"{ani(x,y):>12}" for y in keys))
