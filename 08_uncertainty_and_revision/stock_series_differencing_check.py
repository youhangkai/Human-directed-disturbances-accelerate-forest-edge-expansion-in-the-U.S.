import warnings; warnings.filterwarnings("ignore")
import sys,io,numpy as np,pandas as pd
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding="utf-8")
import statsmodels.api as sm
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.tsa.stattools import adfuller
G=r"G:/Hangkai/CONUS_Forest_Edge_LCMAP"
def st(p): return "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "n.s."

def diff_analysis(name,x,y):
    x=np.asarray(x,float); y=np.asarray(y,float)
    d=np.diff(y); xd=x[1:]; xdc=xd-xd.mean()          # annual increments
    # (a) is the MEAN increment != 0 ?  (AR1 on the increments, intercept only)
    m0=ARIMA(d,order=(1,0,0),trend="c").fit()
    mean_inc,p_mean=m0.params[0],m0.pvalues[0]
    # (b) is there a TREND in the increments ? (= acceleration)
    m1=ARIMA(d,exog=xdc,order=(1,0,0),trend="c").fit()
    acc,p_acc=m1.params[1],m1.pvalues[1]
    # unit-root check on the level series
    try: adf_p=adfuller(y,autolag="AIC")[1]
    except Exception: adf_p=np.nan
    rho_d=np.corrcoef(d[:-1],d[1:])[0,1]
    return dict(series=name,mean_inc=mean_inc,p_mean=p_mean,acc=acc,p_acc=p_acc,
                adf_p=adf_p,rho_diff=rho_d,n=len(d))

rows=[]
a=pd.read_csv(G+r"/key_outputs/forest_area_by_region.csv")
for reg in ["Northeast","Midwest","South","West"]:
    s=a[a.Region==reg].sort_values("Year"); rows.append(diff_analysis(f"Area {reg}",s.Year,s.ForestArea_km2))
ca=a.groupby("Year").ForestArea_km2.sum().reset_index()
rows.append(diff_analysis("Area CONUS",ca.Year,ca.ForestArea_km2))
e=pd.read_csv(G+r"/key_outputs/forest_edge_by_region.csv").sort_values("Year")
col={"Northeast":"Northeast (km)","Midwest":"Midwest (km)","South":"South (km)","West":"West (km)"}
e["CONUS"]=e[[col[r] for r in col]].sum(axis=1)
for reg in ["Northeast","Midwest","South","West"]: rows.append(diff_analysis(f"Edge {reg}",e.Year,e[col[reg]]))
rows.append(diff_analysis("Edge CONUS",e.Year,e.CONUS))
df=pd.DataFrame(rows); df.to_csv("ar1_differenced.csv",index=False)

print("="*116)
print("STOCK VARIABLES RE-TESTED ON ANNUAL INCREMENTS (first differences), AR1 errors")
print("  mean increment != 0  ->  'area declined' / 'edge expanded'")
print("  trend in increments  ->  'the RATE changed' (acceleration)")
print("="*116)
print(f"{'series':<20}{'n':>4}{'ADF p':>8}{'rho(diff)':>11} | {'mean incr/yr':>14}{'p':>10} {'':4}| {'accel/yr2':>11}{'p':>10} {'':4}")
print("-"*116)
for _,r in df.iterrows():
    print(f"{r.series:<20}{int(r.n):>4}{r.adf_p:>8.3f}{r.rho_diff:>11.2f} | {r.mean_inc:>14.4g}{r.p_mean:>10.2e} {st(r.p_mean):<4}| {r.acc:>11.4g}{r.p_acc:>10.2e} {st(r.p_acc):<4}")
print("\nADF p<0.05 => level series is stationary (no unit root). p>0.05 => consistent with a random walk.")
print("saved ar1_differenced.csv")
