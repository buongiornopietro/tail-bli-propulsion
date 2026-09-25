"""Mappa finale: branca turbolenta (solida) e branca a transizione naturale
(fascia tra il criterio di Michel, ottimistico, e il criterio e^N)."""
import en_fig  # noqa: F401  (English figures with FIG_EN=1)
import json, glob, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

def carica(cartella):
    out = {}
    for f in glob.glob(f"{cartella}/soglia_Re*.json"):
        if "seme" in f or f.count("_") > 2: continue
        d = json.load(open(f)); out[(d["Re"], d["regime"])] = d["K_star"]
    return out

eN = carica("out"); mich = carica("out/v3")
Re = sorted({k[0] for k in eN})
lim = 22
fig, ax = plt.subplots(figsize=(9.5, 5.5))
turb = [100 * eN[(r, "forzata")] for r in Re]
ax.plot(Re, turb, "o-", c="#d35400", lw=2.5, ms=7, label="turbolento (superficie ruvida o aria agitata)")
nat_eN = np.array([min(100 * eN[(r, "naturale")], lim) for r in Re])
nat_mi = np.array([min(100 * mich[(r, "naturale")], lim) for r in Re])
ax.fill_between(Re, np.minimum(nat_eN, nat_mi), np.maximum(nat_eN, nat_mi), color="#17a589", alpha=.25,
                label="transizione naturale: fascia di incertezza del criterio")
ax.plot(Re, nat_eN, "s--", c="#17a589", lw=2, ms=6, label="transizione naturale, criterio e^N (stima preferita)")
ax.plot(Re, nat_mi, "^:", c="#148f77", lw=1.5, ms=6, label="transizione naturale, criterio di Michel (ottimistico)")
for r, v in zip(Re, nat_eN):
    if v >= lim: ax.annotate("", xy=(r, lim + 1.5), xytext=(r, lim - 0.5), arrowprops=dict(arrowstyle="->", color="#17a589"))
ax.set_xscale("log"); ax.set_ylim(0, lim + 3); ax.set_xlim(2e5, 1.5e8)
ax.set_xlabel("numero di Reynolds volumetrico  Re_V")
ax.set_ylabel("perdita massima ammissibile nel condotto  K*  [% di ½ρU²]")
ax.set_title("Soglia di convenienza del propulsore in coda\n"
             "sotto la curva conviene il propulsore in coda, sopra l'elica convenzionale", fontsize=10, loc="left")
ax.legend(fontsize=8, loc="upper right"); ax.grid(alpha=.3)
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout(); plt.savefig("out/mappa_finale.png", dpi=120)
print("turbolento:", [f"{v:.1f}" for v in turb])
print("naturale e^N:", [f"{v:.1f}" for v in nat_eN])
print("naturale Michel:", [f"{v:.1f}" for v in nat_mi])
