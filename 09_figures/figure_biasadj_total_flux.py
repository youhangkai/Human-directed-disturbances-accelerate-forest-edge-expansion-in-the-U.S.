import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT=r"C:\Users\hyou34\OneDrive - UW-Madison\manuscripts\paper\CONUS Forest Landscape Dynamics Attribution\Nature version"
# thousand km2
cats=["Gross loss","Gross gain","Net flux\n(gain − loss)"]
mapv =np.array([834.668, 733.623, -101.046])
biav =np.array([447.737, 403.448, -44.289])
biaci=np.array([24.460,  23.552,   33.956])
x=np.arange(3); w=0.38
fig,ax=plt.subplots(figsize=(9.2,5.4))
b1=ax.bar(x-w/2, mapv, w, color="#9aa7b1", label="Map (pixel-counting)")
b2=ax.bar(x+w/2, biav, w, color="#2c6fb3", yerr=biaci, capsize=5,
          error_kw=dict(ecolor="#123", lw=1.3), label="Bias-adjusted (design-based, 95% CI)")
ax.axhline(0,color="k",lw=0.8)
def lab(bars,vals,dy):
    for r,v in zip(bars,vals):
        ax.annotate(f"{v:+.0f}" if v<0 else f"{v:.0f}",
                    (r.get_x()+r.get_width()/2, v),
                    textcoords="offset points", xytext=(0, dy if v>=0 else -dy-6),
                    ha="center", fontsize=9, fontweight="bold")
lab(b1,mapv,6); lab(b2,biav,6)
# cross-check annotation on net
ax.annotate("independent endpoint\npaired estimate: −44.3k\n(agrees)",
            xy=(2+w/2, -44.289), xytext=(1.35, -95),
            fontsize=8.2, style="italic", color="#123", ha="center",
            arrowprops=dict(arrowstyle="->", color="#123", lw=1))
ax.set_xticks(x); ax.set_xticklabels(cats, fontsize=10)
ax.set_ylabel("Forest area flux (thousand km², CONUS 1988–2021)", fontsize=10)
ax.set_title("Total forest flux: map vs. bias-adjusted (design-based), CONUS 1988–2021",
             fontsize=12, fontweight="bold")
ax.legend(fontsize=9.5, frameon=False, loc="lower left")
ax.grid(axis="y", alpha=0.25)
fig.text(0.5,0.03,"Transition-based fluxes summed over 1989–2021 (gross loss = forest→non-forest, gross gain = non-forest→forest).",
         ha="center", fontsize=8.2, style="italic")
fig.text(0.5,0.008,"Gain and net are NOT split by disturbance agent: the reference labels recovery as Growth/Recovery, not by agent.",
         ha="center", fontsize=8.2, style="italic")
plt.tight_layout(rect=[0,0.07,1,1])
out=OUT+r"\Fig_biasadj_total_loss_gain_net.png"
plt.savefig(out,dpi=200); print("saved",out)
