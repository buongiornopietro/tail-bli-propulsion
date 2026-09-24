"""Figura della forma "garantita" (out/robusta.json -> out/robusta.png).
Sinistra: profili delle tre forme per ogni Re_V. Destra: guadagno di velocità a potenza
fissata rispetto alla goccia turbolenta, scenario per scenario."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ris = json.load(open("out/robusta.json"))
COL = {"turbolenta": "0.45", "nominale": "#d35400", "robusta": "#2471a3"}
fig, ax = plt.subplots(len(ris), 2, figsize=(12, 2.6 * len(ris)), squeeze=False,
                       gridspec_kw=dict(width_ratios=[1.3, 1]))
for i, d in enumerate(ris):
    a, b = ax[i]
    for nome, f in d["forme"].items():
        x, r = np.array(f["x"]), np.array(f["r"])
        a.plot(x, r, c=COL[nome], lw=2, label=f"{nome} (L/D {f['LsuD']:.1f})")
        a.plot(x, -r, c=COL[nome], lw=2)
    a.set_aspect("equal"); a.axis("off")
    a.set_title(f"Re_V = {d['ReV']:.0e}", loc="left", fontsize=10)
    a.legend(fontsize=7, loc="lower right", frameon=False)
    sc = list(d["C_rif"])
    xx = np.arange(len(sc)); w = 0.38
    for k, nome in enumerate(("nominale", "robusta")):
        g = [d["forme"][nome]["guadagno_vel"][s] for s in sc]
        b.bar(xx + (k - 0.5) * w, g, w, color=COL[nome], label=nome)
    b.axhline(0, c="0.3", lw=0.8)
    b.set_xticks(xx, sc, fontsize=8); b.set_ylabel("velocità vs goccia [%]", fontsize=8)
    b.grid(axis="y", alpha=.3)
    if i == 0:
        b.legend(fontsize=8, frameon=False)
fig.suptitle("Forma ottima per la transizione naturale: nominale vs garantita (stessa potenza)",
             fontsize=11, x=0.01, ha="left")
plt.tight_layout(); plt.savefig("out/robusta.png", dpi=110)
print("out/robusta.png")
