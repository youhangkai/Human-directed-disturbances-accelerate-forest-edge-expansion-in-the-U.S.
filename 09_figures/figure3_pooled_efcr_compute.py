import sys,io,numpy as np,pandas as pd
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding="utf-8")
from scipy.stats import theilslopes
import pymannkendall as mk
BIG=r"G:\Hangkai\CONUS_Forest_Edge_LCMAP\key_outputs\Forest_Area_Change_historical_Disturbance_Attribution_1988_2021_pa&ownership_reprojected.csv"
OUT=r"C:\Users\hyou34\OneDrive - UW-Madison\manuscripts\paper\CONUS Forest Landscape Dynamics Attribution\Nature version"
REGION={
 'Connecticut':'Northeast','Maine':'Northeast','Massachusetts':'Northeast','New Hampshire':'Northeast','Rhode Island':'Northeast','Vermont':'Northeast','New Jersey':'Northeast','New York':'Northeast','Pennsylvania':'Northeast',
 'Illinois':'Midwest','Indiana':'Midwest','Michigan':'Midwest','Ohio':'Midwest','Wisconsin':'Midwest','Iowa':'Midwest','Kansas':'Midwest','Minnesota':'Midwest','Missouri':'Midwest','Nebraska':'Midwest','North Dakota':'Midwest','South Dakota':'Midwest',
 'Delaware':'South','Florida':'South','Georgia':'South','Maryland':'South','North Carolina':'South','South Carolina':'South','Virginia':'South','West Virginia':'South','District of Columbia':'South','Alabama':'South','Kentucky':'South','Mississippi':'South','Tennessee':'South','Arkansas':'South','Louisiana':'South','Oklahoma':'South','Texas':'South',
 'Arizona':'West','Colorado':'West','Idaho':'West','Montana':'West','Nevada':'West','New Mexico':'West','Utah':'West','Wyoming':'West','California':'West','Oregon':'West','Washington':'West'}
print("reading...",flush=True)
b=pd.read_csv(BIG,usecols=["Ecoregion","Year_From","ForestChangeType","DisturbanceCategory","PixelCount"])
b["Region"]=b.Ecoregion.map(REGION)
b["ct"]=b.ForestChangeType.astype(str).str.strip()
b["fr"]=b.ct.apply(lambda s:'0' if len(s)==1 else s[0]); b["to"]=b.ct.str[-1]
b["frag"]=np.where((b.fr=='5')&(b.to.isin(['1','2','3','4'])),b.PixelCount,0)
b["loss"]=np.where((b.fr.isin(['1','2','3','4','5']))&(b.to=='0'),b.PixelCount,0)
b["D"]=b.DisturbanceCategory.replace({"Forest Management":"Logging","Other":"Others"})
b=b[(b.Year_From>=1988)&(b.Year_From<=2020)]
AG=["Stress","Logging","Others","Water Dynamic","Natural Hazard","Agriculture Activity","Construction","Fire"]

def stats(sub):
    """pooled EFCR per year within this spatial unit, then mean + Sen slope + MK p"""
    g=sub.groupby("Year_From")[["frag","loss"]].sum()
    g=g[g.loss>0]
    if len(g)<8: return None
    y=(g.frag/g.loss).values; x=g.index.values.astype(float)
    sl,ic,lo,hi=theilslopes(y,x,0.95); r=mk.original_test(y)
    return dict(mean=y.mean(),sen=sl,lo=lo,hi=hi,p=r.p,n=len(y),
                loss_km2=g.loss.sum()*0.0009)

# ---- panel (a): national pooled ----
rows=[]
for d in AG:
    s=stats(b[b.D==d])
    if s: rows.append(dict(Agent=d,**s))
A=pd.DataFrame(rows).sort_values("mean",ascending=False)
A.to_csv(OUT+r"\Fig3a_pooled_EFCR_by_disturbance.csv",index=False)
print("\n=== PANEL (a) national POOLED EFCR ===")
print(A[["Agent","mean","sen","lo","hi","p","loss_km2"]].round(4).to_string(index=False))

# ---- panel (b): region x disturbance ----
rows=[]
for d in AG:
    for r_ in ["Northeast","Midwest","South","West"]:
        s=stats(b[(b.D==d)&(b.Region==r_)])
        if s: rows.append(dict(Agent=d,Region=r_,**s))
B=pd.DataFrame(rows)
B.to_csv(OUT+r"\Fig3b_pooled_EFCR_region_disturbance.csv",index=False)
print(f"\n=== PANEL (b) {len(B)} region x disturbance points ===")
print(B[["Agent","Region","mean","sen","p","loss_km2"]].round(4).to_string(index=False))
print("\nsaved both CSVs")
