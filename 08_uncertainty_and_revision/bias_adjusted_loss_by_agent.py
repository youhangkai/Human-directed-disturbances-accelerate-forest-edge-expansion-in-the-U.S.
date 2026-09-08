import sys,io,numpy as np,pandas as pd
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding="utf-8")
REF=r"G:\Hangkai\CONUS_Forest_Edge_LCMAP\bias_adjusted_area\ref\LCMAP_Collection1.3_simple.xlsx"
BIG=r"G:\Hangkai\CONUS_Forest_Edge_LCMAP\key_outputs\Forest_Area_Change_historical_Disturbance_Attribution_1988_2021_pa&ownership_reprojected.csv"
OUT=r"C:\Users\hyou34\OneDrive - UW-Madison\manuscripts\paper\CONUS Forest Landscape Dynamics Attribution\Nature version"
SCR=r"C:\Users\hyou34\AppData\Local\Temp\claude\C--Users-hyou34\1af3a46b-9cd5-47d3-849e-d949261f3d6b\scratchpad"
A=8081894.0; Z=1.959964; KM2=0.0009
AGENTS=["Logging","Fire","Stress","Construction+Agriculture","Water Dynamic","Natural Hazard","Others"]

# ---------- 1) reference loss transitions ----------
print("reading ref xlsx...",flush=True)
df=pd.read_excel(REF,usecols=["plotid","image_year","LCMAP","change_process"]).sort_values(["plotid","image_year"])
n=df.plotid.nunique(); gg=df.groupby("plotid")
df["plc"]=gg["LCMAP"].shift(1); df["pyr"]=gg["image_year"].shift(1); df["pcp"]=gg["change_process"].shift(1)
loss=df[(df.plc=="Tree Cover")&(df.LCMAP!="Tree Cover")&(df.pyr==df.image_year-1)&(df.image_year>=1989)&(df.image_year<=2021)].copy()
def ag(r):
    a=r.change_process
    if pd.isna(a) or a in("Stable","Growth/Recovery"):
        return r.pcp if (pd.notna(r.pcp) and r.pcp not in("Stable","Growth/Recovery")) else a
    return a
loss["agent"]=loss.apply(ag,axis=1)
xw={"Harvest":"Logging","Fire":"Fire","Mechanical":"Construction+Agriculture","Structural Decline":"Stress","Spectral Decline":"Stress","Hydrology":"Water Dynamic","Wind":"Natural Hazard","Debris":"Natural Hazard"}
loss["G"]=loss.agent.map(lambda a: xw.get(a,"Others"))
loss[["plotid","image_year","agent","G"]].to_csv(SCR+r"\loss_transitions.csv",index=False)
print(f"ref loss transitions: {len(loss)} (n_plots={n}); saved",flush=True)

# ---------- 2) map loss by agent from big csv ----------
print("reading 2.9GB attribution csv...",flush=True)
b=pd.read_csv(BIG,usecols=["Year_From","ForestChangeType","DisturbanceCategory","PixelCount"])
b["fr"]=b.ForestChangeType.astype(str).str[0]; b["to"]=b.ForestChangeType.astype(str).str[1]
b=b[(b.fr.isin(['1','2','3','4','5']))&(b.to=='0')]   # forest loss
b["D"]=b.DisturbanceCategory.replace({"No Disturbance":"Others","Forest Management":"Logging","Other":"Others","No Disturbance Detected":"Others"})
mmap={"Logging":"Logging","Fire":"Fire","Stress":"Stress","Construction":"Construction+Agriculture","Agriculture Activity":"Construction+Agriculture","Water Dynamic":"Water Dynamic","Natural Hazard":"Natural Hazard"}
b["G"]=b.D.map(lambda d: mmap.get(d,"Others"))
mapyr=b.groupby(["Year_From","G"])["PixelCount"].sum().reset_index()
mapyr["map_km2"]=mapyr.PixelCount*KM2
mapyr.to_csv(SCR+r"\map_loss_by_agent_annual.csv",index=False)
def mapc(g): return mapyr[mapyr.G==g].map_km2.sum()
def mapy(g,yf):
    s=mapyr[(mapyr.G==g)&(mapyr.Year_From==yf)].map_km2.sum(); return s

# ---------- 3) cumulative table ----------
rows=[]
for g in AGENTS:
    sub=loss[loss.G==g]; ne=len(sub); cpp=sub.groupby("plotid").size(); full=np.zeros(n); full[:len(cpp)]=cpp.values
    rows.append(dict(Agent=g,n_events=ne,Map_loss_km2=round(mapc(g)),BiasAdj_loss_km2=round(ne/n*A),CI95_km2=round(Z*A*np.sqrt(full.var(ddof=1)/n))))
cpp=loss.groupby("plotid").size(); full=np.zeros(n); full[:len(cpp)]=cpp.values
rows.append(dict(Agent="TOTAL",n_events=len(loss),Map_loss_km2=round(sum(mapc(g) for g in AGENTS)),BiasAdj_loss_km2=round(len(loss)/n*A),CI95_km2=round(Z*A*np.sqrt(full.var(ddof=1)/n))))
cum=pd.DataFrame(rows); cum.to_csv(OUT+r"\TableSx_loss_by_agent_cumulative.csv",index=False)
print("\n===== CUMULATIVE 1988-2021 =====\n"+cum.to_string(index=False))

# ---------- 4) annual table ----------
arr=[]
for yf in range(1988,2021):
    ty=yf+1
    for g in AGENTS:
        ne=int(((loss.G==g)&(loss.image_year==ty)).sum()); p=ne/n
        arr.append(dict(Year_From=yf,Agent=g,Map_loss_km2=round(mapy(g,yf)),BiasAdj_loss_km2=round(p*A),CI95_km2=round(Z*A*np.sqrt(p*(1-p)/n)),n_events=ne))
ann=pd.DataFrame(arr); ann.to_csv(OUT+r"\TableSx_loss_by_agent_annual.csv",index=False)
print("\n===== ANNUAL saved. Logging & Fire preview (BiasAdj km2) =====")
print(ann[ann.Agent.isin(["Logging","Fire"])].pivot(index="Year_From",columns="Agent",values="BiasAdj_loss_km2").to_string())
print("\nSAVED both to manuscript folder.")
