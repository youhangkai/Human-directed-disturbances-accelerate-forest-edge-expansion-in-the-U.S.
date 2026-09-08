import sys,io,numpy as np,pandas as pd
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding="utf-8")
BIG=r"G:\Hangkai\CONUS_Forest_Edge_LCMAP\key_outputs\Forest_Area_Change_historical_Disturbance_Attribution_1988_2021_pa&ownership_reprojected.csv"
OUT=r"C:\Users\hyou34\OneDrive - UW-Madison\manuscripts\paper\CONUS Forest Landscape Dynamics Attribution\Nature version"
KM2=0.0009
REGION={
 'Connecticut':'Northeast','Maine':'Northeast','Massachusetts':'Northeast','New Hampshire':'Northeast','Rhode Island':'Northeast','Vermont':'Northeast','New Jersey':'Northeast','New York':'Northeast','Pennsylvania':'Northeast',
 'Illinois':'Midwest','Indiana':'Midwest','Michigan':'Midwest','Ohio':'Midwest','Wisconsin':'Midwest','Iowa':'Midwest','Kansas':'Midwest','Minnesota':'Midwest','Missouri':'Midwest','Nebraska':'Midwest','North Dakota':'Midwest','South Dakota':'Midwest',
 'Delaware':'South','Florida':'South','Georgia':'South','Maryland':'South','North Carolina':'South','South Carolina':'South','Virginia':'South','West Virginia':'South','District of Columbia':'South','Alabama':'South','Kentucky':'South','Mississippi':'South','Tennessee':'South','Arkansas':'South','Louisiana':'South','Oklahoma':'South','Texas':'South',
 'Arizona':'West','Colorado':'West','Idaho':'West','Montana':'West','Nevada':'West','New Mexico':'West','Utah':'West','Wyoming':'West','California':'West','Oregon':'West','Washington':'West'}
print("reading big csv...",flush=True)
b=pd.read_csv(BIG,usecols=["Ecoregion","Year_From","ForestChangeType","DisturbanceCategory","PixelCount"])
b["Region"]=b.Ecoregion.map(REGION)
print("unmapped states:",b[b.Region.isna()].Ecoregion.unique())
b["ct"]=b.ForestChangeType.astype(str).str.strip()
b["fr"]=b.ct.apply(lambda s:'0' if len(s)==1 else s[0]); b["to"]=b.ct.str[-1]
mmap={"Forest Management":"Logging","Fire":"Fire","Stress":"Stress","Construction":"Construction+Agriculture",
      "Agriculture Activity":"Construction+Agriculture","Water Dynamic":"Water Dynamic","Natural Hazard":"Natural Hazard"}
b["G"]=b.DisturbanceCategory.map(lambda d:mmap.get(d,"Others"))
b["loss_px"]=np.where((b.fr.isin(list('12345')))&(b.to=='0'),b.PixelCount,0)
b["gain_px"]=np.where((b.fr=='0')&(b.to.isin(list('12345'))),b.PixelCount,0)

# regional totals (loss/gain/net, all agents)
rt=b.groupby("Region").agg(loss=("loss_px","sum"),gain=("gain_px","sum")).reset_index()
rt["loss_km2"]=rt.loss*KM2; rt["gain_km2"]=rt.gain*KM2; rt["net_km2"]=rt.gain_km2-rt.loss_km2
rt=rt[["Region","loss_km2","gain_km2","net_km2"]]
conus=pd.DataFrame([dict(Region="CONUS",loss_km2=rt.loss_km2.sum(),gain_km2=rt.gain_km2.sum(),net_km2=rt.net_km2.sum())])
rt=pd.concat([rt,conus],ignore_index=True)
rt.to_csv(OUT+r"\Table_region_loss_gain_net_map.csv",index=False)
print("\n=== REGIONAL map loss/gain/net (km2) ===\n"+rt.round().to_string(index=False))

# regional loss by agent (map) + dominant agent
ra=b.groupby(["Region","G"]).agg(loss_km2=("loss_px",lambda s:s.sum()*KM2)).reset_index()
ra.to_csv(OUT+r"\Table_region_loss_by_agent_map.csv",index=False)
print("\n=== dominant loss agent by region (map) ===")
for reg in ["Northeast","Midwest","South","West"]:
    s=ra[ra.Region==reg].sort_values("loss_km2",ascending=False)
    tot=s.loss_km2.sum()
    top=s.iloc[0]
    print(f"{reg:10s}: total loss {tot:8.0f} km2 | top = {top.G} {top.loss_km2:.0f} ({100*top.loss_km2/tot:.0f}%)")
print("\nSAVED regional tables.")
