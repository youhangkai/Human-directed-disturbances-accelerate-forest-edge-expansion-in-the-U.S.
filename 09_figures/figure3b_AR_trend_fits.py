import warnings; warnings.filterwarnings("ignore")
import sys,io,numpy as np,pandas as pd
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding="utf-8")
import statsmodels.api as sm
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.stats.diagnostic import acorr_ljungbox
NV=r"C:/Users/hyou34/OneDrive - UW-Madison/manuscripts/paper/CONUS Forest Landscape Dynamics Attribution/Nature version"
SUB=NV+"/Statistical_Revision_AR"

def fit(x,y,maxp=3):
    x=np.asarray(x,float); y=np.asarray(y,float)
    ok=np.isfinite(y); x=x[ok]; y=y[ok]; n=len(y)
    if n<10: return dict(n=n,slope=np.nan,se=np.nan,lo=np.nan,hi=np.nan,p=np.nan,ar=np.nan,lb=np.nan,mean=np.nanmean(y) if n else np.nan)
    tc=x-x.mean(); best=None
    for p in range(0,maxp+1):
        try:
            m=ARIMA(y,exog=tc,order=(p,0,0),trend="c").fit()
            if best is None or m.aic<best[1].aic: best=(p,m)
        except Exception: pass
    p,m=best; ci=m.conf_int()[1]
    try: lb=acorr_ljungbox(m.resid,lags=[5],return_df=True)["lb_pvalue"].values[0]
    except Exception: lb=np.nan
    return dict(n=n,slope=m.params[1],se=m.bse[1],lo=ci[0],hi=ci[1],p=m.pvalues[1],ar=p,lb=lb,mean=y.mean())

rows=[]
# 1) region x disturbance (for Fig 3b)
g=pd.read_csv(SUB+"/efcr_annual_region_x_disturbance.csv")
for (r,d),s in g.groupby(["Region","D"]):
    s=s.sort_values("Year_From"); f=fit(s.Year_From,s.EFCR)
    rows.append(dict(level="region x disturbance",Region=r,Disturbance=d,**f))
# 2) national by disturbance
nd=pd.read_csv(NV+"/efcr_annual_by_disturbance.csv")
for d,s in nd.groupby("D"):
    s=s.sort_values("Year_From"); f=fit(s.Year_From,s.EFCR)
    rows.append(dict(level="national by disturbance",Region="CONUS",Disturbance=d,**f))
# 3) regional aggregate + CONUS, full period and 2001-2020
for r in ["Northeast","Midwest","South","West","CONUS"]:
    s=pd.read_csv(NV+f"/efcr_annual_{r}.csv").sort_values("Year_From"); s["EFCR"]=s.frag/s.loss
    rows.append(dict(level="regional aggregate 1988-2020",Region=r,Disturbance="All",**fit(s.Year_From,s.EFCR)))
    s2=s[s.Year_From>=2001]
    rows.append(dict(level="regional aggregate 2001-2020",Region=r,Disturbance="All",**fit(s2.Year_From,s2.EFCR,maxp=1)))
df=pd.DataFrame(rows)
df["sig"]=df.p<0.05
df.to_csv(SUB+"/AR_trend_results_all_series.csv",index=False)
print("rows:",len(df))
print("\nregion x disturbance fits with n<10 (skipped):")
print(df[(df.level=="region x disturbance")&(df.n<10)][["Region","Disturbance","n"]].to_string(index=False) if (df[(df.level=="region x disturbance")&(df.n<10)]).shape[0] else "  none")
print("\nLjung-Box p<0.05 (residual dependence not removed):")
bad=df[(df.lb<0.05)]
print(bad[["level","Region","Disturbance","ar","lb"]].to_string(index=False) if len(bad) else "  none")
print("\nsaved AR_trend_results_all_series.csv")
