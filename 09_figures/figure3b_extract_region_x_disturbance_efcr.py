import warnings; warnings.filterwarnings("ignore")
import sys,io,numpy as np,pandas as pd
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding="utf-8")
BIG=r"G:/Hangkai/CONUS_Forest_Edge_LCMAP/key_outputs/Forest_Area_Change_historical_Disturbance_Attribution_1988_2021_pa&ownership_reprojected.csv"
SUB=r"C:/Users/hyou34/OneDrive - UW-Madison/manuscripts/paper/CONUS Forest Landscape Dynamics Attribution/Nature version/Statistical_Revision_AR"
REGION={**{s:'Northeast' for s in ['Connecticut','Maine','Massachusetts','New Hampshire','Rhode Island','Vermont','New Jersey','New York','Pennsylvania']},
 **{s:'Midwest' for s in ['Illinois','Indiana','Michigan','Ohio','Wisconsin','Iowa','Kansas','Minnesota','Missouri','Nebraska','North Dakota','South Dakota']},
 **{s:'South' for s in ['Delaware','Florida','Georgia','Maryland','North Carolina','South Carolina','Virginia','West Virginia','District of Columbia','Alabama','Kentucky','Mississippi','Tennessee','Arkansas','Louisiana','Oklahoma','Texas']},
 **{s:'West' for s in ['Arizona','Colorado','Idaho','Montana','Nevada','New Mexico','Utah','Wyoming','California','Oregon','Washington']}}
print("reading big csv...",flush=True)
b=pd.read_csv(BIG,usecols=["Ecoregion","Year_From","ForestChangeType","DisturbanceCategory","PixelCount"])
b["Region"]=b.Ecoregion.map(REGION)
b["ct"]=b.ForestChangeType.astype(str).str.strip()
b["fr"]=b.ct.apply(lambda s:'0' if len(s)==1 else s[0]); b["to"]=b.ct.str[-1]
b["frag"]=np.where((b.fr=='5')&(b.to.isin(['1','2','3','4'])),b.PixelCount,0)
b["loss"]=np.where((b.fr.isin(['1','2','3','4','5']))&(b.to=='0'),b.PixelCount,0)
b["D"]=b.DisturbanceCategory.replace({"Forest Management":"Logging","Other":"Others"})
b=b[(b.Year_From>=1988)&(b.Year_From<=2020)]
g=b.groupby(["Region","D","Year_From"])[["frag","loss"]].sum().reset_index()
g["EFCR"]=np.where(g.loss>0,g.frag/g.loss,np.nan)
g.to_csv(SUB+"/efcr_annual_region_x_disturbance.csv",index=False)
print("saved region x disturbance annual series:",g.Region.nunique(),"regions x",g.D.nunique(),"agents")
