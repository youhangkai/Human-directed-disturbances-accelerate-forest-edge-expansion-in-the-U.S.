import warnings; warnings.filterwarnings("ignore")
import sys,io,numpy as np,pandas as pd
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding="utf-8")
from statsmodels.tsa.arima.model import ARIMA
NV=r"C:/Users/hyou34/OneDrive - UW-Madison/manuscripts/paper/CONUS Forest Landscape Dynamics Attribution/Nature version"
SUB=NV+"/Statistical_Revision_AR"

def sel(y,x,maxp,crit="aic"):
    tc=x-x.mean(); best=None
    for p in range(0,maxp+1):
        try:
            m=ARIMA(y,exog=tc,order=(p,0,0),trend="c").fit()
            k=len(m.params); n=len(y)
            val=m.aic if crit=="aic" else m.aic+2*k*(k+1)/max(n-k-1,1)
            if best is None or val<best[0]: best=(val,p,m)
        except Exception: pass
    _,p,m=best; return p,m.params[1],m.pvalues[1],m.conf_int()[1]

def sig(p): return "sig" if p<0.05 else "n.s."
print("="*104)
print("SENSITIVITY OF CONCLUSIONS TO THE AR-ORDER CEILING")
print("="*104)

# --- EFCR regional, full and recent ---
print("\n[A] Regional aggregate EFCR")
print(f"{'series':<26}{'ceiling':>8}{'AIC p*':>7}{'slope':>10}{'pval':>10}{'':5}{'AICc p*':>8}{'slope':>10}{'pval':>10}")
for r in ["Northeast","Midwest","South","West","CONUS"]:
    d=pd.read_csv(NV+f"/efcr_annual_{r}.csv").sort_values("Year_From"); d["E"]=d.frag/d.loss
    for lab,sub in [("1988-2020",d),("2001-2020",d[d.Year_From>=2001])]:
        y=sub.E.values.astype(float); x=sub.Year_From.values.astype(float)
        row=f"{r+' '+lab:<26}"
        for c in ([1,3] if lab=="2001-2020" else [3,5]):
            p,b,pv,ci=sel(y,x,c); row+=f"{c:>8}{p:>7}{b:>10.4f}{pv:>10.3f} {sig(pv):<4}"
        pa,ba,pva,cia=sel(y,x,3 if lab=="2001-2020" else 5,crit="aicc")
        row+=f"{pa:>8}{ba:>10.4f}{pva:>10.3f}"
        print(row)

# --- acceleration ---
print("\n[B] Edge acceleration gamma")
e=pd.read_csv(r"G:/Hangkai/CONUS_Forest_Edge_LCMAP/key_outputs/forest_edge_by_region.csv").sort_values("Year")
col={r:f"{r} (km)" for r in ["Northeast","Midwest","South","West"]}
e["CONUS"]=e[[col[r] for r in col]].sum(axis=1)
print(f"{'series':<18}{'ceil=3':>22}{'ceil=4':>22}{'ceil=5':>22}")
for r in ["Northeast","Midwest","South","West","CONUS"]:
    yv=(e[col[r]] if r!="CONUS" else e["CONUS"]).values.astype(float)
    dd=np.diff(yv); xd=e.Year.values[1:].astype(float)
    row=f"{'Edge '+r:<18}"
    for c in [3,4,5]:
        p,b,pv,ci=sel(dd,xd,c); row+=f"  AR({p}) {b:>8.0f} p={pv:<8.1e}"
    print(row)

# --- region x disturbance: do any markers flip? ---
print("\n[C] Figure 3b markers: ceiling 3 vs ceiling 5")
g=pd.read_csv(SUB+"/efcr_annual_region_x_disturbance.csv")
flips=[]
for (rg,dz),s in g.groupby(["Region","D"]):
    if dz=="No Disturbance": continue
    s=s.sort_values("Year_From"); s=s[s.loss>0]
    y=s.EFCR.values.astype(float); x=s.Year_From.values.astype(float)
    p3,b3,pv3,_=sel(y,x,3); p5,b5,pv5,_=sel(y,x,5)
    if (pv3<=0.05)!=(pv5<=0.05): flips.append((rg,dz,p3,pv3,p5,pv5))
print(f"  markers changing significance: {len(flips)} of 32")
for rg,dz,p3,pv3,p5,pv5 in flips: print(f"    {rg} x {dz}: AR({p3}) p={pv3:.3f} -> AR({p5}) p={pv5:.3f}")
