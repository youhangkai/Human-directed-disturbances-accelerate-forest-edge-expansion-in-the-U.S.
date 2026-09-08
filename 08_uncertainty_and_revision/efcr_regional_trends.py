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
print("reading big csv...",flush=True)
b=pd.read_csv(BIG,usecols=["Ecoregion","Year_From","ForestChangeType","PixelCount"])
b["Region"]=b.Ecoregion.map(REGION)
b["ct"]=b.ForestChangeType.astype(str).str.strip()
b["fr"]=b.ct.apply(lambda s:'0' if len(s)==1 else s[0]); b["to"]=b.ct.str[-1]
b["frag_px"]=np.where((b.fr=='5')&(b.to.isin(['1','2','3','4'])),b.PixelCount,0)   # interior->exterior
b["loss_px"]=np.where((b.fr.isin(['1','2','3','4','5']))&(b.to=='0'),b.PixelCount,0) # forest loss

def series(sub):
    g=sub.groupby("Year_From").agg(frag=("frag_px","sum"),loss=("loss_px","sum")).reset_index()
    g=g[(g.Year_From>=1988)&(g.Year_From<=2020)].sort_values("Year_From")   # 33 yrs, matches prior n=33
    g["EFCR"]=g.frag/g.loss
    return g

def trend(g):
    x=g.Year_From.values.astype(float); y=g.EFCR.values.astype(float)
    sl,ic,lo,hi=theilslopes(y,x,0.95)
    r=mk.original_test(y)
    return sl,lo,hi,r.p,r.Tau,y[0],y[-1],np.nanmean(y)

rows=[]
for reg in ["Northeast","Midwest","South","West"]:
    g=series(b[b.Region==reg]); sl,lo,hi,p,tau,e0,e1,mn=trend(g)
    rows.append(dict(Region=reg,EFCR_1988=e0,EFCR_2020=e1,EFCR_mean=mn,Sen_per_yr=sl,CI_lo=lo,CI_hi=hi,MK_p=p,Tau=tau))
    g.to_csv(OUT+rf"\efcr_annual_{reg}.csv",index=False)
gC=series(b); sl,lo,hi,p,tau,e0,e1,mn=trend(gC)
rows.append(dict(Region="CONUS",EFCR_1988=e0,EFCR_2020=e1,EFCR_mean=mn,Sen_per_yr=sl,CI_lo=lo,CI_hi=hi,MK_p=p,Tau=tau))
gC.to_csv(OUT+r"\efcr_annual_CONUS.csv",index=False)
res=pd.DataFrame(rows)
res.to_csv(OUT+r"\Table_EFCR_trends_by_region.csv",index=False)
pd.set_option("display.width",200)
sh=res.copy()
for c in ["EFCR_1988","EFCR_2020","EFCR_mean","Sen_per_yr","CI_lo","CI_hi","Tau"]:
    sh[c]=sh[c].map(lambda v:f"{v:.4f}")
sh["MK_p"]=res.MK_p.map(lambda v:f"{v:.2e}")
sh["sig"]=res.MK_p.map(lambda v:"***" if v<0.001 else "**" if v<0.01 else "*" if v<0.05 else "n.s.")
print(sh.to_string(index=False))
print("\nsaved Table_EFCR_trends_by_region.csv + per-region annual series")
