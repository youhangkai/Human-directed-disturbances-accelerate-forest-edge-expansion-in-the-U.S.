import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT=r"C:\Users\hyou34\OneDrive - UW-Madison\manuscripts\paper\CONUS Forest Landscape Dynamics Attribution\Nature version"
d=pd.read_csv(OUT+r"\Table_loss_gain_net_by_agent_annual.csv")
d["yr"]=d.Year_From+1
order=[("Logging","Logging"),("Fire","Fire"),("Construction+Agriculture","Construction/Agriculture"),
       ("Stress","Stress"),("Water Dynamic","Water Dynamic"),("Natural Hazard","Natural Hazard"),("Others","Others")]
fig,axes=plt.subplots(2,4,figsize=(15,6.4),sharex=True); axes=axes.ravel()
for k,(ag,lab) in enumerate(order):
    ax=axes[k]; s=d[d.G==ag].sort_values("yr")
    x=s.yr.values; loss=s.loss_km2.values/1e3; gain=s.gain_km2.values/1e3; net=s.net_km2.values/1e3
    ax.bar(x,-loss,color="#c44",width=0.9,label="Gross loss")
    ax.bar(x, gain,color="#4a9",width=0.9,label="Gross gain")
    ax.plot(x,net,"-",color="k",lw=1.4,label="Net flux")
    ax.axhline(0,color="k",lw=0.6)
    ax.set_title(lab, fontsize=10, fontweight="bold")
    ax.grid(alpha=0.22); ax.margins(x=0)
    if k%4==0: ax.set_ylabel("Forest flux\n(thousand km²/yr)",fontsize=8.5)
axes[7].axis("off")
axes[7].legend(*axes[0].get_legend_handles_labels(),loc="center",fontsize=11,frameon=False)
for ax in axes[4:7]: ax.set_xlabel("Year",fontsize=8.5)
fig.suptitle("Annual forest loss (down), gain (up), and net flux by disturbance agent, CONUS 1989–2021 (map pixel-counting)",
             fontsize=11.5, fontweight="bold", y=0.995)
fig.text(0.5,0.005,"Loss = forest→non-forest, gain = non-forest→forest, attributed to the pixel's mapped disturbance history. Net flux = gain − loss.",
         ha="center",fontsize=8.5,style="italic")
plt.tight_layout(rect=[0,0.02,1,0.965])
out=OUT+r"\Fig_annual_loss_gain_net.png"
plt.savefig(out,dpi=180); print("saved",out)
