import sys,io,numpy as np,pandas as pd
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding="utf-8")
BIG=r"G:\Hangkai\CONUS_Forest_Edge_LCMAP\key_outputs\Forest_Area_Change_historical_Disturbance_Attribution_1988_2021_pa&ownership_reprojected.csv"
OUT=r"C:\Users\hyou34\OneDrive - UW-Madison\manuscripts\paper\CONUS Forest Landscape Dynamics Attribution\Nature version"
KM2=0.0009
print("reading big csv...",flush=True)
b=pd.read_csv(BIG,usecols=["Year_From","ForestChangeType","DisturbanceCategory","PixelCount"])
b["ct"]=b.ForestChangeType.astype(str).str.strip()
b["fr"]=b.ct.apply(lambda s:'0' if len(s)==1 else s[0]); b["to"]=b.ct.str[-1]
mmap={"Forest Management":"Logging","Fire":"Fire","Stress":"Stress","Construction":"Construction+Agriculture",
      "Agriculture Activity":"Construction+Agriculture","Water Dynamic":"Water Dynamic","Natural Hazard":"Natural Hazard"}
b["G"]=b.DisturbanceCategory.map(lambda d:mmap.get(d,"Others"))
b["loss_px"]=np.where((b.fr.isin(list('12345')))&(b.to=='0'),b.PixelCount,0)
b["gain_px"]=np.where((b.fr=='0')&(b.to.isin(list('12345'))),b.PixelCount,0)
g=b.groupby(["Year_From","G"]).agg(loss=("loss_px","sum"),gain=("gain_px","sum")).reset_index()
g["loss_km2"]=g.loss*KM2; g["gain_km2"]=g.gain*KM2; g["net_km2"]=g.gain_km2-g.loss_km2
g.to_csv(OUT+r"\Table_loss_gain_net_by_agent_annual.csv",index=False)
cum=g.groupby("G").agg(loss_km2=("loss_km2","sum"),gain_km2=("gain_km2","sum")).reset_index()
cum["net_km2"]=cum.gain_km2-cum.loss_km2
cum.to_csv(OUT+r"\Table_loss_gain_net_by_agent_cumulative.csv",index=False)
print("\n=== CUMULATIVE (km2) ===\n"+cum.round().to_string(index=False))
print("\nTOTAL loss",round(cum.loss_km2.sum()),"gain",round(cum.gain_km2.sum()),"net",round(cum.net_km2.sum()))
print("saved annual + cumulative CSVs")
