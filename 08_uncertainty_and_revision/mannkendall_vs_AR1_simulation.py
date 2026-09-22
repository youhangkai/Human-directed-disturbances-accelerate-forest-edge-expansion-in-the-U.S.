import warnings; warnings.filterwarnings("ignore")
import sys,io,numpy as np,pandas as pd
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding="utf-8")
from scipy.stats import theilslopes
import pymannkendall as mk
import statsmodels.api as sm
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.stats.stattools import durbin_watson

G=r"G:/Hangkai/CONUS_Forest_Edge_LCMAP"
OUT=r"C:/Users/hyou34/OneDrive - UW-Madison/manuscripts/paper/CONUS Forest Landscape Dynamics Attribution/Nature version"
def stars(p): return "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "n.s."

def analyse(name,x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); n=len(y); xc=x-x.mean()
    r=mk.original_test(y); sen,_,_,_=theilslopes(y,x,0.95)
    X=sm.add_constant(xc); ols=sm.OLS(y,X).fit(); resid=ols.resid
    rho=np.corrcoef(resid[:-1],resid[1:])[0,1]; dw=durbin_watson(resid)
    m=ARIMA(y,exog=xc,order=(1,0,0),trend="c").fit()
    # params order: [const, x1, ar.L1, sigma2]
    b1,se,p_ar1,phi=m.params[1],m.bse[1],m.pvalues[1],m.params[2]
    phi_lo,phi_hi=m.conf_int()[2]
    return dict(series=name,n=n,sen=sen,p_mk=r.p,p_ols=ols.pvalues[1],rho=rho,dw=dw,
                ar1_slope=b1,ar1_se=se,p_ar1=p_ar1,phi=phi,phi_lo=phi_lo,phi_hi=phi_hi)

rows=[]
a=pd.read_csv(G+r"/key_outputs/forest_area_by_region.csv")
for reg in ["Northeast","Midwest","South","West"]:
    s=a[a.Region==reg].sort_values("Year"); rows.append(analyse(f"Area {reg}",s.Year,s.ForestArea_km2))
ca=a.groupby("Year").ForestArea_km2.sum().reset_index()
rows.append(analyse("Area CONUS",ca.Year,ca.ForestArea_km2))
b=pd.read_csv(G+r"/bias_adjusted_area/bias_adjusted_forest_area.csv").dropna(subset=["adj_area_km2"]).sort_values("year")
rows.append(analyse("Area CONUS bias-adj",b.year,b.adj_area_km2))
e=pd.read_csv(G+r"/key_outputs/forest_edge_by_region.csv").sort_values("Year")
col={"Northeast":"Northeast (km)","Midwest":"Midwest (km)","South":"South (km)","West":"West (km)"}
e["CONUS"]=e[[col[r] for r in col]].sum(axis=1)
for reg in ["Northeast","Midwest","South","West"]: rows.append(analyse(f"Edge {reg}",e.Year,e[col[reg]]))
rows.append(analyse("Edge CONUS",e.Year,e.CONUS))
for reg in ["Northeast","Midwest","South","West","CONUS"]:
    d=pd.read_csv(OUT+rf"/efcr_annual_{reg}.csv").sort_values("Year_From")
    rows.append(analyse(f"EFCR {reg}",d.Year_From,d.frag/d.loss))
df=pd.DataFrame(rows); df.to_csv("ar1_vs_mk.csv",index=False)

pd.set_option("display.width",250)
print("="*122)
print("MANN-KENDALL vs AR1 REGRESSION      rho = lag-1 autocorr of OLS resid | phi = fitted AR1 coef")
print("="*122)
print(f"{'series':<22}{'n':>4}{'rho':>7}{'DW':>6} | {'MK p':>10} {'':4}| {'AR1 slope':>12}{'± SE':>11} {'AR1 p':>10} {'':4}| {'phi [95% CI]':>20}  change")
print("-"*122)
for _,r in df.iterrows():
    chg=""
    if r.p_mk<0.05 and r.p_ar1>=0.05: chg="LOST sig."
    elif r.p_mk>=0.05 and r.p_ar1<0.05: chg="GAINED sig."
    print(f"{r.series:<22}{int(r.n):>4}{r.rho:>7.2f}{r.dw:>6.2f} | {r.p_mk:>10.2e} {stars(r.p_mk):<4}| {r.ar1_slope:>12.4g}{r.ar1_se:>11.4g} {r.p_ar1:>10.2e} {stars(r.p_ar1):<4}| {r.phi:>6.2f} [{r.phi_lo:5.2f},{r.phi_hi:5.2f}]  {chg}")

# ---------- false-positive simulation calibrated to OUR autocorrelation ----------
print("\n"+"="*122)
print("FALSE-POSITIVE SIMULATION: zero-trend AR1 noise, n=34, 5000 reps. How often does each test claim p<0.05?")
print("="*122)
rng=np.random.default_rng(42); N=34; R=5000
print(f"{'phi':>6}{'MK false +':>14}{'OLS false +':>14}{'AR1 false +':>14}   (nominal 5%)")
print("-"*60)
for phi in [0.0,0.3,0.5,0.7,0.9,0.95]:
    fm=fo=fa=0
    for _ in range(R):
        w=rng.standard_normal(N); y=np.zeros(N)
        for t in range(1,N): y[t]=phi*y[t-1]+w[t]
        x=np.arange(N,dtype=float); xc=x-x.mean()
        if mk.original_test(y).p<0.05: fm+=1
        if sm.OLS(y,sm.add_constant(xc)).fit().pvalues[1]<0.05: fo+=1
        try:
            mm=ARIMA(y,exog=xc,order=(1,0,0),trend="c").fit()
            if mm.pvalues[1]<0.05: fa+=1
        except Exception: pass
    print(f"{phi:>6.2f}{100*fm/R:>13.1f}%{100*fo/R:>13.1f}%{100*fa/R:>13.1f}%")
print("\nsaved ar1_vs_mk.csv")
