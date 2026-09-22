import warnings; warnings.filterwarnings("ignore")
import sys,io,numpy as np,pandas as pd
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding="utf-8")
import statsmodels.api as sm
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.stats.diagnostic import acorr_ljungbox
G=r"G:/Hangkai/CONUS_Forest_Edge_LCMAP"

def acf1(r): return np.corrcoef(r[:-1],r[1:])[0,1]

def accel(name,x,y):
    """ dE_t = alpha + gamma*t + eps_t ;  choose AR(p) for eps by AIC; also Newey-West HAC """
    x=np.asarray(x,float); y=np.asarray(y,float)
    d=np.diff(y); t=x[1:]; tc=t-t.mean(); n=len(d)
    X=sm.add_constant(tc)
    ols=sm.OLS(d,X).fit()
    rho_d=acf1(ols.resid)
    # --- AR(p) selection by AIC ---
    best=None
    for p in range(0,5):
        try:
            m=ARIMA(d,exog=tc,order=(p,0,0),trend="c").fit()
            if best is None or m.aic<best[1].aic: best=(p,m)
        except Exception: pass
    p,m=best
    g=m.params[1]; se=m.bse[1]; pv=m.pvalues[1]
    ci=m.conf_int()[1]
    # residual adequacy of chosen model
    lb=acorr_ljungbox(m.resid,lags=[5],return_df=True)["lb_pvalue"].values[0]
    rho_after=acf1(m.resid)
    # --- Newey-West HAC (conservative, no AR structure assumed) ---
    L=int(np.floor(4*(n/100)**(2/9)))
    hac=sm.OLS(d,X).fit(cov_type="HAC",cov_kwds={"maxlags":max(L,1)})
    g_h,se_h,p_h=hac.params[1],hac.bse[1],hac.pvalues[1]
    ci_h=hac.conf_int()[1]
    # --- cross-check: quadratic on levels, gamma should = 2c ---
    tl=x-x.mean(); Xq=np.column_stack([np.ones_like(tl),tl,tl**2])
    q=sm.OLS(y,Xq).fit(); two_c=2*q.params[2]
    return dict(series=name,n=n,rho_diff=rho_d,ar_p=p,gamma=g,se=se,lo=ci[0],hi=ci[1],p=pv,
                rho_after=rho_after,lb_p=lb,g_hac=g_h,se_hac=se_h,lo_h=ci_h[0],hi_h=ci_h[1],p_hac=p_h,
                hac_L=max(L,1),two_c=two_c)

rows=[]
e=pd.read_csv(G+r"/key_outputs/forest_edge_by_region.csv").sort_values("Year")
col={"Northeast":"Northeast (km)","Midwest":"Midwest (km)","South":"South (km)","West":"West (km)"}
e["CONUS"]=e[[col[r] for r in col]].sum(axis=1)
for reg in ["Northeast","Midwest","South","West"]: rows.append(accel(f"Edge {reg}",e.Year,e[col[reg]]))
rows.append(accel("Edge CONUS",e.Year,e.CONUS))
a=pd.read_csv(G+r"/key_outputs/forest_area_by_region.csv")
ca=a.groupby("Year").ForestArea_km2.sum().reset_index()
rows.append(accel("Area CONUS",ca.Year,ca.ForestArea_km2))
df=pd.DataFrame(rows); df.to_csv("accel_proper.csv",index=False)

def st(p): return "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "n.s."
print("="*126)
print("ACCELERATION  gamma  from   dE_t = alpha + gamma*t + eps_t     (gamma = annual change in the rate of edge expansion)")
print("="*126)
print(f"{'series':<15}{'n':>3}{'rho(d)':>8}{'AR':>4} | {'gamma':>10}{'SE':>9} {'95% CI':>20} {'p':>9} {'':4}| {'LB p':>6}{'rho_res':>8} | {'gamma HAC':>10}{'p HAC':>9} {'':4}| {'2c quad':>10}")
print("-"*126)
for _,r in df.iterrows():
    print(f"{r.series:<15}{int(r.n):>3}{r.rho_diff:>8.2f}{int(r.ar_p):>4} | {r.gamma:>10.0f}{r.se:>9.0f} [{r.lo:>8.0f},{r.hi:>8.0f}] {r.p:>9.1e} {st(r.p):<4}| {r.lb_p:>6.2f}{r.rho_after:>8.2f} | {r.g_hac:>10.0f}{r.p_hac:>9.1e} {st(r.p_hac):<4}| {r.two_c:>10.0f}")
print(f"\nHAC maxlag used = {int(df.hac_L.iloc[0])};  LB p = Ljung-Box(5) on chosen-model residuals (want > 0.05)")
print("'2c quad' = 2 x quadratic coefficient from levels; should match gamma if the two parameterisations agree.")
