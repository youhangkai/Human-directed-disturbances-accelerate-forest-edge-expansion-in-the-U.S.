import sys,io,numpy as np,pandas as pd
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding="utf-8")
from scipy.stats import theilslopes
try:
    import pymannkendall as mk
    HAVEMK=True
except Exception:
    HAVEMK=False
G=r"G:\Hangkai\CONUS_Forest_Edge_LCMAP"
OUT=r"C:\Users\hyou34\OneDrive - UW-Madison\manuscripts\paper\CONUS Forest Landscape Dynamics Attribution\Nature version"

area=pd.read_csv(G+r"\key_outputs\forest_area_by_region.csv")           # Year,Region,ForestArea_km2
edge=pd.read_csv(G+r"\key_outputs\forest_edge_by_region.csv")           # Year, <Region> (km)
ba=pd.read_csv(G+r"\bias_adjusted_area\bias_adjusted_forest_area.csv")

def trend(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float)
    sl,ic,lo,hi=theilslopes(y,x,0.95)
    if HAVEMK:
        r=mk.original_test(y); p=r.p; tau=r.Tau
    else:
        p=np.nan; tau=np.nan
    return sl,lo,hi,p,tau

REGS=["Northeast","Midwest","South","West"]
rows=[]
# ---- AREA (map) per region ----
for reg in REGS:
    s=area[area.Region==reg].sort_values("Year")
    y1=s[s.Year==1988].ForestArea_km2.values[0]; y2=s[s.Year==2021].ForestArea_km2.values[0]
    sl,lo,hi,p,tau=trend(s.Year,s.ForestArea_km2)
    rows.append(dict(Metric="Forest area (map, km²)",Region=reg,V1988=y1,V2021=y2,Change=y2-y1,
                     PctChange=100*(y2-y1)/y1,Sen_per_yr=sl,Sen_lo=lo,Sen_hi=hi,MK_p=p))
# CONUS area map (sum of regions)
ca=area.groupby("Year").ForestArea_km2.sum().reset_index()
y1=ca[ca.Year==1988].ForestArea_km2.values[0]; y2=ca[ca.Year==2021].ForestArea_km2.values[0]
sl,lo,hi,p,tau=trend(ca.Year,ca.ForestArea_km2)
rows.append(dict(Metric="Forest area (map, km²)",Region="CONUS",V1988=y1,V2021=y2,Change=y2-y1,
                 PctChange=100*(y2-y1)/y1,Sen_per_yr=sl,Sen_lo=lo,Sen_hi=hi,MK_p=p))

# ---- AREA (bias-adjusted, CONUS only, post-stratified) ----
bb=ba.dropna(subset=["adj_area_km2"]).sort_values("year")
y1=bb[bb.year==1988].adj_area_km2.values[0]; y2=bb[bb.year==2021].adj_area_km2.values[0]
c1=bb[bb.year==1988].adj_ci95_km2.values[0]; c2=bb[bb.year==2021].adj_ci95_km2.values[0]
sl,lo,hi,p,tau=trend(bb.year,bb.adj_area_km2)
rows.append(dict(Metric="Forest area (bias-adj, km²)",Region="CONUS",V1988=y1,V2021=y2,Change=y2-y1,
                 PctChange=100*(y2-y1)/y1,Sen_per_yr=sl,Sen_lo=lo,Sen_hi=hi,MK_p=p,
                 CI1988=c1,CI2021=c2))

# ---- EDGE length (map) per region ----
ecol={"Northeast":"Northeast (km)","Midwest":"Midwest (km)","South":"South (km)","West":"West (km)"}
edge["CONUS (km)"]=edge[[ecol[r] for r in REGS]].sum(axis=1)
for reg in REGS+["CONUS"]:
    col=ecol.get(reg,"CONUS (km)")
    s=edge.sort_values("Year"); yy=s[col].values; xx=s.Year.values
    y1=s[s.Year==1988][col].values[0]; y2=s[s.Year==2021][col].values[0]
    sl,lo,hi,p,tau=trend(xx,yy)
    rows.append(dict(Metric="Forest edge length (map, km)",Region=reg,V1988=y1,V2021=y2,Change=y2-y1,
                     PctChange=100*(y2-y1)/y1,Sen_per_yr=sl,Sen_lo=lo,Sen_hi=hi,MK_p=p))

df=pd.DataFrame(rows)
df.to_csv(OUT+r"\Table_area_edge_trends.csv",index=False)
pd.set_option("display.width",200,"display.max_columns",20)
def fmt(v):
    return f"{v:,.0f}" if abs(v)>=100 else f"{v:,.2f}"
show=df.copy()
for c in ["V1988","V2021","Change","Sen_per_yr","Sen_lo","Sen_hi"]:
    show[c]=show[c].map(lambda v: fmt(v) if pd.notna(v) else "")
show["PctChange"]=show.PctChange.map(lambda v:f"{v:+.1f}%" if pd.notna(v) else "")
show["MK_p"]=show.MK_p.map(lambda v:f"{v:.1e}" if pd.notna(v) else "")
print(show[["Metric","Region","V1988","V2021","Change","PctChange","Sen_per_yr","Sen_lo","Sen_hi","MK_p"]].to_string(index=False))
print("\nsaved Table_area_edge_trends.csv")
