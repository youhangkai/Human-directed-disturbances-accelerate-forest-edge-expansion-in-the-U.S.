import pandas as pd, numpy as np, sys, io
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding="utf-8")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT=r"C:\Users\hyou34\OneDrive - UW-Madison\manuscripts\paper\CONUS Forest Landscape Dynamics Attribution\Nature version"
A=pd.read_csv(OUT+r"\Fig3a_pooled_EFCR_by_disturbance.csv")     # CONUS by agent
B=pd.read_csv(OUT+r"\Fig3b_pooled_EFCR_region_disturbance.csv") # region x agent
R=pd.read_csv(OUT+r"\Table_EFCR_trends_by_region.csv")          # all-disturbance by region + CONUS

REGS=["Northeast","Midwest","South","West"]
AG=list(A.sort_values("mean",ascending=False).Agent)            # pooled CONUS order

# ---- assemble matrix: rows = agent (+ All), cols = regions + CONUS ----
rows=[]
for ag in AG:
    d={"Disturbance":ag}
    for r in REGS:
        m=B[(B.Agent==ag)&(B.Region==r)]
        d[r]=m["mean"].values[0] if len(m) else np.nan
        d[r+"_p"]=m["p"].values[0] if len(m) else np.nan
    d["CONUS"]=A[A.Agent==ag]["mean"].values[0]
    d["CONUS_p"]=A[A.Agent==ag]["p"].values[0]
    rows.append(d)
d={"Disturbance":"All disturbances"}
for r in REGS:
    d[r]=R[R.Region==r].EFCR_mean.values[0]; d[r+"_p"]=R[R.Region==r].MK_p.values[0]
d["CONUS"]=R[R.Region=="CONUS"].EFCR_mean.values[0]; d["CONUS_p"]=R[R.Region=="CONUS"].MK_p.values[0]
rows.append(d)
M=pd.DataFrame(rows)
M.to_csv(OUT+r"\Table_EFCR_mean_by_region_and_CONUS.csv",index=False)

print("="*94)
print("MEAN POOLED EFCR BY REGION AND CONUS (1988-2020, n=33)   * = MK trend p<=0.05")
print("="*94)
hdr=f"{'Disturbance':<22}"+"".join(f"{r:>13}" for r in REGS)+f"{'CONUS':>13}"
print(hdr); print("-"*94)
for _,r in M.iterrows():
    line=f"{r.Disturbance:<22}"
    for c in REGS+["CONUS"]:
        star="*" if r[c+"_p"]<=0.05 else " "
        line+=f"{r[c]:>12.2f}{star}"
    print(line)

# ---------------- FIGURE ----------------
COLR={"Northeast":"#1b7837","Midwest":"#7fbc41","South":"#d6604d","West":"#4393c3"}
lab=[a.replace(" Activity","").replace("Natural Hazard","Natural\nHazard").replace("Water Dynamic","Water\nDynamic")
     .replace("All disturbances","ALL\ndisturbances") for a in M.Disturbance]
n=len(M); x=np.arange(n); w=0.16
fig,ax=plt.subplots(figsize=(15,6.2))
for i,r in enumerate(REGS):
    ax.bar(x+(i-2)*w, M[r], w, color=COLR[r], label=r, edgecolor="white", linewidth=0.6)
ax.bar(x+2*w, M["CONUS"], w, color="#222222", label="CONUS (pooled)", edgecolor="white", linewidth=0.6)
for xi,v in zip(x+2*w, M["CONUS"]):
    ax.annotate(f"{v:.2f}",(xi,v),textcoords="offset points",xytext=(0,3),ha="center",fontsize=7.6,fontweight="bold")
ax.axvline(n-1.5,color="k",lw=0.9,ls="--",alpha=0.5)
ax.set_xticks(x); ax.set_xticklabels(lab,fontsize=9.2)
ax.set_ylabel("Mean EFCR  (interior→exterior area per unit forest loss)",fontsize=10.5)
ax.set_title("Mean EFCR by disturbance type: four Census regions vs. CONUS (area-weighted, 1988–2020)",
             fontsize=12,fontweight="bold")
ax.legend(fontsize=9.5,frameon=False,ncol=5,loc="upper right")
ax.grid(axis="y",alpha=0.25); ax.set_axisbelow(True)
ax.set_ylim(0,max(M[REGS].max().max(),M.CONUS.max())*1.15)
fig.text(0.5,0.008,"CONUS bar is the area-weighted (pooled) value, not the average of the four regional bars: regions with little disturbance area contribute little to it. "
         "Bars right of the dashed line pool all disturbance types.",ha="center",fontsize=8.4,style="italic")
plt.tight_layout(rect=[0,0.035,1,1])
out=OUT+r"\Fig_EFCR_by_region_and_CONUS.png"
plt.savefig(out,dpi=220); print("\nsaved",out)
print("saved Table_EFCR_mean_by_region_and_CONUS.csv")
