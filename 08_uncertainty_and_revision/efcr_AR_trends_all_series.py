import warnings; warnings.filterwarnings("ignore")
import sys,io,numpy as np,pandas as pd
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding="utf-8")
import statsmodels.api as sm
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.stats.diagnostic import acorr_ljungbox
import pymannkendall as mk
from scipy.stats import theilslopes
OUT=r"C:/Users/hyou34/OneDrive - UW-Madison/manuscripts/paper/CONUS Forest Landscape Dynamics Attribution/Nature version"
def st(p): return "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "n.s."

def fit(name,x,y,maxp=None):
    x=np.asarray(x,float); y=np.asarray(y,float); n=len(y); tc=x-x.mean()
    if n<8: return None
    if maxp is None: maxp=1 if n<15 else 3
    mkr=mk.original_test(y); sen,_,_,_=theilslopes(y,x,0.95)
    ols=sm.OLS(y,sm.add_constant(tc)).fit(); rho=np.corrcoef(ols.resid[:-1],ols.resid[1:])[0,1]
    best=None
    for p in range(0,maxp+1):
        try:
            m=ARIMA(y,exog=tc,order=(p,0,0),trend="c").fit()
            if best is None or m.aic<best[1].aic: best=(p,m)
        except Exception: pass
    if best is None: return None
    p,m=best; b=m.params[1]; se=m.bse[1]; pv=m.pvalues[1]; ci=m.conf_int()[1]
    try: lb=acorr_ljungbox(m.resid,lags=[min(5,n//3)],return_df=True)["lb_pvalue"].values[0]
    except Exception: lb=np.nan
    rho_a=np.corrcoef(m.resid[:-1],m.resid[1:])[0,1]
    return dict(series=name,n=n,rho=rho,ar_p=p,sen=sen,p_mk=mkr.p,slope=b,se=se,lo=ci[0],hi=ci[1],
                p_ar=pv,lb=lb,rho_after=rho_a)

rows=[]
# ---- regional aggregate EFCR, full + subperiods ----
for reg in ["Northeast","Midwest","South","West","CONUS"]:
    d=pd.read_csv(OUT+rf"/efcr_annual_{reg}.csv").sort_values("Year_From")
    d["E"]=d.frag/d.loss
    for lab,sub in [("1988-2020",d),("2001-2020",d[d.Year_From>=2001]),("2011-2020",d[d.Year_From>=2011])]:
        r=fit(f"{reg} {lab}",sub.Year_From,sub.E)
        if r: rows.append(r)
# ---- national EFCR by disturbance type ----
import os
f=OUT+"/efcr_annual_by_disturbance.csv"
if os.path.exists(f):
    g=pd.read_csv(f)
    for ag,sub in g.groupby("D"):
        sub=sub.sort_values("Year_From")
        r=fit(f"[{ag}] 1988-2020",sub.Year_From,sub.EFCR)
        if r: rows.append(r)
else:
    print("!! efcr_annual_by_disturbance.csv not ready yet\n")
df=pd.DataFrame(rows); df.to_csv("efcr_ar1_full.csv",index=False)
print("="*128)
print("EFCR TRENDS: Mann-Kendall vs AR(p) regression, with residual diagnostics")
print("="*128)
print(f"{'series':<26}{'n':>4}{'rho':>7}{'AR':>4} | {'Sen':>9}{'MK p':>10} {'':4}| {'AR slope':>10}{'SE':>9} {'95% CI':>20}{'AR p':>10} {'':4}| {'LB':>5}{'rho_a':>7}")
print("-"*128)
prev=""
for _,r in df.iterrows():
    tag=r.series.split()[0]
    if tag!=prev and prev: print("-"*128)
    prev=tag
    flag=""
    if r.p_mk<0.05 and r.p_ar>=0.05: flag=" LOST"
    elif r.p_mk>=0.05 and r.p_ar<0.05: flag=" GAIN"
    warn="!" if (not np.isnan(r.lb) and r.lb<0.05) else " "
    print(f"{r.series:<26}{int(r.n):>4}{r.rho:>7.2f}{int(r.ar_p):>4} | {r.sen:>9.4f}{r.p_mk:>10.1e} {st(r.p_mk):<4}| {r.slope:>10.4f}{r.se:>9.4f} [{r.lo:>8.4f},{r.hi:>8.4f}]{r.p_ar:>10.1e} {st(r.p_ar):<4}| {r.lb:>5.2f}{r.rho_after:>7.2f}{warn}{flag}")
print("\n'!' = Ljung-Box p<0.05, residual dependence NOT adequately removed -> inference unreliable")
