import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
OUT=r"C:\Users\hyou34\OneDrive - UW-Madison\manuscripts\paper\CONUS Forest Landscape Dynamics Attribution\Nature version"
A=pd.read_csv(OUT+r"\Fig3a_pooled_EFCR_by_disturbance.csv")
B=pd.read_csv(OUT+r"\Fig3b_pooled_EFCR_region_disturbance.csv")

# --- original palettes, kept per panel as in the published figure ---
# panel (a) bar palette
COL_A={"Stress":"#F49AC1","Logging":"#4FC3A1","Natural Hazard":"#A9603F","Others":"#B3B3B3",
       "Water Dynamic":"#4A90D9","Agriculture Activity":"#F2C14E","Construction":"#7E3F98","Fire":"#E63329"}
# panel (b) scatter palette (matches the original legend)
COL={"Logging":"#3FA34D","Construction":"#6B2D8B","Stress":"#E63329","Natural Hazard":"#F2D024",
     "Water Dynamic":"#2E86D4","Fire":"#F06AAE","Agriculture Activity":"#F5A623","Others":"#9E9E9E"}
MRK={"Midwest":"o","Northeast":"s","South":"^","West":"D"}
SHORT={"Agriculture Activity":"Agriculture\nActivity","Natural Hazard":"Natural\nHazard",
       "Water Dynamic":"Water\nDynamic","Construction":"Construction","Logging":"Logging",
       "Stress":"Stress","Fire":"Fire","Others":"Others"}

fig=plt.figure(figsize=(14,5.6))
gs=fig.add_gridspec(1,2,width_ratios=[1,1.25],wspace=0.22)

# ---------------- panel (a) ----------------
ax=fig.add_subplot(gs[0,0])
A=A.sort_values("mean",ascending=False)
x=np.arange(len(A))
ax.bar(x,A["mean"],color=[COL_A[a] for a in A.Agent],edgecolor="white",linewidth=0.8)
for xi,v in zip(x,A["mean"]):
    ax.annotate(f"{v:.2f}",(xi,v),textcoords="offset points",xytext=(0,3),
                ha="center",fontsize=8.5,fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels([SHORT[a] for a in A.Agent],rotation=40,ha="right",fontsize=9)
ax.set_ylabel("Mean EFCR",fontsize=10.5)
ax.set_title("(a) Mean EFCR by disturbance type",fontsize=11,fontweight="bold",loc="left")
ax.grid(axis="y",alpha=0.25); ax.set_axisbelow(True)
ax.set_ylim(0,max(A["mean"])*1.18)

# ---------------- panel (b) ----------------
ax2=fig.add_subplot(gs[0,1])
for _,r in B.iterrows():
    sig=r.p<=0.05
    ax2.scatter(r["mean"],r.sen,marker=MRK[r.Region],s=78,
                facecolor=COL[r.Agent] if sig else "none",
                edgecolor=COL[r.Agent],linewidth=1.5,zorder=3)
ax2.axhline(0,color="k",lw=0.7,zorder=1)
ax2.set_xlabel("Mean EFCR",fontsize=10.5)
ax2.set_ylabel("Sen's slope of EFCR (yr$^{-1}$)",fontsize=10.5)
ax2.set_title("(b) EFCR and its trend by disturbance type and region",fontsize=11,fontweight="bold",loc="left")
ax2.grid(alpha=0.25); ax2.set_axisbelow(True)

h_reg=[Line2D([],[],marker=MRK[k],color="k",ls="none",ms=8,mfc="none",label=k) for k in ["Midwest","Northeast","South","West"]]
h_dis=[Line2D([],[],marker="o",color=COL[k],ls="none",ms=8,label=k) for k in
       ["Logging","Construction","Stress","Natural Hazard","Water Dynamic","Fire","Agriculture Activity","Others"]]
h_sig=[Line2D([],[],marker="o",color="k",ls="none",ms=8,mfc="k",label="Significant (p ≤ 0.05)"),
       Line2D([],[],marker="o",color="k",ls="none",ms=8,mfc="none",label="Not significant")]
l1=ax2.legend(handles=h_reg,title="Region",loc="upper left",bbox_to_anchor=(1.02,1.0),fontsize=8.5,title_fontsize=9,frameon=False)
l2=ax2.legend(handles=h_dis,title="Disturbance",loc="upper left",bbox_to_anchor=(1.02,0.72),fontsize=8.5,title_fontsize=9,frameon=False)
l3=ax2.legend(handles=h_sig,title="Trend significance",loc="upper left",bbox_to_anchor=(1.02,0.24),fontsize=8.5,title_fontsize=9,frameon=False)
ax2.add_artist(l1); ax2.add_artist(l2)

fig.suptitle("Figure 3. Exterior Forest Conversion Ratio (EFCR) by disturbance type and region — area-weighted (pooled)",
             fontsize=12,fontweight="bold",y=0.995)
plt.tight_layout(rect=[0,0,0.86,0.95])
out=OUT+r"\Figure_3_pooled.png"
plt.savefig(out,dpi=250,bbox_inches="tight"); print("saved",out)
