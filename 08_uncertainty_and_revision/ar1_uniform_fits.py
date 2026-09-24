import warnings; warnings.filterwarnings("ignore")
import sys,io,numpy as np,pandas as pd
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding="utf-8")
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.stats.diagnostic import acorr_ljungbox
NV=r"C:/Users/hyou34/OneDrive - UW-Madison/manuscripts/paper/CONUS Forest Landscape Dynamics Attribution/Nature version"
SUB=NV+"/Statistical_Revision_AR"; G=r"G:/Hangkai/CONUS_Forest_Edge_LCMAP"

def ar1(x,y):
    """regression on year with AR(1) errors; returns coef, se, CI, p, phi, diagnostics"""
    x=np.asarray(x,float); y=np.asarray(y,float)
    ok=np.isfinite(y); x,y=x[ok],y[ok]; n=len(y)
    if n<10: return None
    tc=x-x.mean()
    m=ARIMA(y,exog=tc,order=(1,0,0),trend="c").fit()
    b,se,p=m.params[1],m.bse[1],m.pvalues[1]; ci=m.conf_int()[1]; phi=m.params[2]
    lag=min(5,max(1,n//3))
    try: lb=acorr_ljungbox(m.resid,lags=[lag],return_df=True)["lb_pvalue"].values[0]
    except Exception: lb=np.nan
    rho=np.corrcoef(m.resid[:-1],m.resid[1:])[0,1]
    return dict(n=n,mean=y.mean(),slope=b,se=se,lo=ci[0],hi=ci[1],p=p,phi=phi,lb=lb,resid_rho=rho)

rows=[]
# 1) region x disturbance (Figure 3b)
g=pd.read_csv(SUB+"/efcr_annual_region_x_disturbance.csv")
for (r,d),s in g.groupby(["Region","D"]):
    s=s.sort_values("Year_From"); f=ar1(s.Year_From,s.EFCR)
    if f: rows.append(dict(level="region x disturbance",Region=r,Disturbance=d,**f))
# 2) national by disturbance
nd=pd.read_csv(NV+"/efcr_annual_by_disturbance.csv")
for d,s in nd.groupby("D"):
    s=s.sort_values("Year_From"); f=ar1(s.Year_From,s.EFCR)
    if f: rows.append(dict(level="national by disturbance",Region="CONUS",Disturbance=d,**f))
# 3) regional aggregate, full + recent
for r in ["Northeast","Midwest","South","West","CONUS"]:
    s=pd.read_csv(NV+f"/efcr_annual_{r}.csv").sort_values("Year_From"); s["E"]=s.frag/s.loss
    for lab,sub in [("EFCR regional 1988-2021",s),("EFCR regional 2001-2021",s[s.Year_From>=2001])]:
        f=ar1(sub.Year_From,sub.E)
        if f: rows.append(dict(level=lab,Region=r,Disturbance="All",**f))
# 4) acceleration on annual increments
e=pd.read_csv(G+"/key_outputs/forest_edge_by_region.csv").sort_values("Year")
col={r:f"{r} (km)" for r in ["Northeast","Midwest","South","West"]}
e["CONUS"]=e[[col[r] for r in col]].sum(axis=1)
for r in ["Northeast","Midwest","South","West","CONUS"]:
    yv=(e[col[r]] if r!="CONUS" else e["CONUS"]).values.astype(float)
    f=ar1(e.Year.values[1:],np.diff(yv))
    if f: rows.append(dict(level="Edge acceleration (km yr-2)",Region=r,Disturbance="Edge",**f))
a=pd.read_csv(G+"/key_outputs/forest_area_by_region.csv").groupby("Year").ForestArea_km2.sum().reset_index()
f=ar1(a.Year.values[1:],np.diff(a.ForestArea_km2.values.astype(float)))
if f: rows.append(dict(level="Area acceleration (km2 yr-2)",Region="CONUS",Disturbance="Area",**f))

df=pd.DataFrame(rows); df["sig"]=df.p<0.05
df.to_csv(SUB+"/AR1_uniform_results.csv",index=False)
print("fits:",len(df))
bad=df[df.lb<0.05]
print("\nLjung-Box p<0.05 (AR(1) did NOT fully remove residual dependence):",len(bad))
if len(bad): print(bad[["level","Region","Disturbance","n","phi","lb","resid_rho"]].to_string(index=False))
print("\nphi range: %.2f to %.2f"%(df.phi.min(),df.phi.max()))
