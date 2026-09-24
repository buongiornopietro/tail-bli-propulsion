import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
B = json.load(open("out/bli_risultati.json")); A = json.load(open("out/axi_risultati.json"))
fig = plt.figure(figsize=(13, 9)); gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.1])
for col, sn, tit in [(0, "A_turbolento_ReV1e7", "Re_V = 10^7, turbolento"), (1, "C_naturale_ReV1e6", "Re_V = 10^6, transizione naturale")]:
    ax = fig.add_subplot(gs[0, col])
    o = A[sn]["ottimo"]; x, r = np.array(o["x"]), np.array(o["r"])
    ax.plot(x, r, c="#7f8c8d", lw=1.5, label=f"ottimo con elica (L/D {o['f']:.1f})"); ax.plot(x, -r, c="#7f8c8d", lw=1.5)
    for vn, c in [("ideale", "#2471a3"), ("realistica", "#d35400")]:
        v = B[f"{sn}|{vn}"]; x, r = np.array(v["x"]), np.array(v["r"])
        ax.plot(x, r, c=c, lw=2.2, label=f"con propulsore in coda, {vn} (L/D {v['f']:.1f})"); ax.plot(x, -r, c=c, lw=2.2)
        xs = v["xs"] * x[-1]; rs = np.interp(xs, x, r)
        ax.plot([xs, xs], [rs, rs + 0.06], c=c, lw=3); ax.plot([xs, xs], [-rs, -rs - 0.06], c=c, lw=3)
    ax.set_aspect("equal"); ax.legend(fontsize=7.5, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    ax.set_title(tit + "\n(trattini = fessura di aspirazione)", fontsize=10, loc="left", pad=62)
    ax.spines[["top", "right"]].set_visible(False)
    ax = fig.add_subplot(gs[1, col])
    labs, vals, cols = [], [], []
    for vn in ("ideale", "realistica"):
        v = B[f"{sn}|{vn}"]
        best_bli = min(v["CP"], v["CP_bli_su_corpo_conv"])
        labs += [f"elica classica\n({vn})", f"propulsore in coda\n({vn})"]
        vals += [v["rif_convenzionale"], best_bli]; cols += ["#7f8c8d", "#2471a3" if vn == "ideale" else "#d35400"]
    b = ax.bar(labs, vals, color=cols)
    for i in (1, 3):
        ax.text(i, vals[i], f"{100*(vals[i]/vals[i-1]-1):+.1f}%", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.set_ylabel("coefficiente di potenza C_P (più basso = meglio)"); ax.tick_params(axis="x", labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("Propulsore in coda (ingestione dello strato limite) vs elica classica — potenza a pari volume e velocità", fontsize=12)
plt.tight_layout(); plt.savefig("out/propulsore_coda.png", dpi=110); print("ok")
