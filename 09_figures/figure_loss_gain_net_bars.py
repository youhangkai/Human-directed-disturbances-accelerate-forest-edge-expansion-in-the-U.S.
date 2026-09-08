import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT=r"C:\Users\hyou34\OneDrive - UW-Madison\manuscripts\paper\CONUS Forest Landscape Dynamics Attribution\Nature version"
c=pd.read_csv(OUT+r"\Table_loss_gain_net_by_agent_cumulative.csv")
label={"Logging":"Logging","Fire":"Fire","Construction+Agriculture":"Construction/\nAgriculture",
       "Stress":"Stress","Water Dynamic":"Water\nDynamic","Natural Hazard":"Natural\nHazard","Others":"Others"}
order=["Logging","Fire","Construction+Agriculture","Stress","Water Dynamic","Natural Hazard","Others"]
c=c.set_index("G").loc[order].reset_index()
x=np.arange(len(order)); w=0.38
fig,ax=plt.subplots(figsize=(11,5.2))
ax.bar(x-w/2, c.loss_km2/1e3, w, color="#c44",label="Gross loss")
ax.bar(x+w/2, c.gain_km2/1e3, w, color="#4a9",label="Gross gain")
ax.plot(x, c.net_km2/1e3, "ko-", ms=6, lw=1.6, label="Net flux (gain − loss)")
for xi,nv in zip(x,c.net_km2/1e3):
    ax.annotate(f"{nv:+.0f}", (xi, nv), textcoords="offset points", xytext=(0,-14 if nv<0 else 8),
                ha="center", fontsize=8.5, fontweight="bold")
ax.axhline(0,color="k",lw=0.7)
ax.set_xticks(x); ax.set_xticklabels([label[o] for o in order], fontsize=9.5)
ax.set_ylabel("Forest area flux (thousand km², 1988–2021)", fontsize=10)
ax.legend(fontsize=9.5, frameon=False, loc="upper right")
ax.set_title("Forest loss, gain, and net flux by disturbance agent, CONUS 1988–2021 (map pixel-counting)",
             fontsize=11.5, fontweight="bold")
ax.grid(axis="y", alpha=0.25)
fig.text(0.5,0.01,"Gross loss = forest→non-forest, gross gain = non-forest→forest, each attributed to the pixel's mapped disturbance history. "
         "Net flux = gain − loss (transition events, not endpoint area change).",
         ha="center", fontsize=8.2, style="italic")
plt.tight_layout(rect=[0,0.035,1,1])
out=OUT+r"\Fig_loss_gain_net_cumulative.png"
plt.savefig(out,dpi=200); print("saved",out)
