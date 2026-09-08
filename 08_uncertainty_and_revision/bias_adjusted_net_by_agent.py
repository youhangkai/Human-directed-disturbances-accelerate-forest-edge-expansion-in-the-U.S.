import sys,io,numpy as np,pandas as pd
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding="utf-8")
REF=r"G:\Hangkai\CONUS_Forest_Edge_LCMAP\bias_adjusted_area\ref\LCMAP_Collection1.3_simple.xlsx"
OUT=r"C:\Users\hyou34\OneDrive - UW-Madison\manuscripts\paper\CONUS Forest Landscape Dynamics Attribution\Nature version"
A=8081894.0; Z=1.959964
AGENTS=["Logging","Fire","Construction+Agriculture","Stress","Water Dynamic","Natural Hazard","Others"]
xw={"Harvest":"Logging","Fire":"Fire","Mechanical":"Construction+Agriculture","Structural Decline":"Stress",
    "Spectral Decline":"Stress","Hydrology":"Water Dynamic","Wind":"Natural Hazard","Debris":"Natural Hazard"}
def grp(a): return xw.get(a,"Others")

print("reading ref xlsx...",flush=True)
df=pd.read_excel(REF,usecols=["plotid","image_year","LCMAP","change_process"]).sort_values(["plotid","image_year"])
N=df.plotid.nunique(); gg=df.groupby("plotid")
df["plc"]=gg["LCMAP"].shift(1); df["pyr"]=gg["image_year"].shift(1); df["pcp"]=gg["change_process"].shift(1)
val=(df.pyr==df.image_year-1)&(df.image_year>=1989)&(df.image_year<=2021)

# ---- LOSS events: agent from change_process (fallback prev-year), same as loss table ----
lm=val&(df.plc=="Tree Cover")&(df.LCMAP!="Tree Cover")
loss=df[lm].copy()
def lag(r):
    a=r.change_process
    if pd.isna(a) or a in("Stable","Growth/Recovery"):
        return r.pcp if (pd.notna(r.pcp) and r.pcp not in("Stable","Growth/Recovery")) else a
    return a
loss["G"]=loss.apply(lambda r: grp(lag(r)),axis=1)
loss["yr"]=loss.image_year.astype(int)

# ---- build per-plot ordered loss-event list for tracing ----
loss_by_plot={}
for pid,sub in loss.groupby("plotid"):
    loss_by_plot[pid]=list(zip(sub.yr.tolist(),sub.G.tolist()))

# ---- GAIN events: attribute to MOST RECENT PRIOR loss agent on same plot (map-analog history) ----
gm=val&(df.plc!="Tree Cover")&(df.LCMAP=="Tree Cover")
gain=df[gm].copy(); gain["yr"]=gain.image_year.astype(int)
def trace(r):
    evs=loss_by_plot.get(r.plotid,[])
    prior=[g for (y,g) in evs if y<r.yr]
    return prior[-1] if prior else "Others"
gain["G"]=gain.apply(trace,axis=1)

# ---- per-plot signed net counts: +1 gain, -1 loss, per (plot,agent,year) ----
le=loss[["plotid","yr","G"]].assign(s=-1)
ge=gain[["plotid","yr","G"]].assign(s=+1)
ev=pd.concat([le,ge],ignore_index=True)

def est_series(agent):
    rows=[]
    sub=ev[ev.G==agent]
    for y in range(1989,2022):
        e=sub[sub.yr==y]
        # per-plot net for this agent-year
        perplot=e.groupby("plotid").s.sum()
        s=perplot.sum(); ss=(perplot**2).sum(); nz=len(perplot)
        mean=s/N
        # variance over all N plots (zeros for plots with no event)
        var=(ss - s*s/N)/(N-1)
        area=mean*A; ci=Z*A*np.sqrt(var/N)
        # also gross loss/gain for reference
        gl=(e.s<0).sum(); gp=(e.s>0).sum()
        rows.append(dict(yr=y,Agent=agent,net_km2=area,ci_km2=ci,n_loss=int(gl),n_gain=int(gp)))
    return rows

out=[]
for g in AGENTS: out+=est_series(g)
res=pd.DataFrame(out)
res.to_csv(OUT+r"\Table_biasadj_net_by_agent_annual.csv",index=False)

print("\n=== cumulative check (bias-adjusted net km2 by agent) ===")
cum=res.groupby("Agent").apply(lambda d:pd.Series(dict(
    net=d.net_km2.sum(), n_loss=d.n_loss.sum(), n_gain=d.n_gain.sum()))).reset_index()
print(cum.round().to_string(index=False))
print("TOTAL net:",round(res.net_km2.sum()))
print("saved Table_biasadj_net_by_agent_annual.csv")
