import en_fig  # noqa: F401  (English figures with FIG_EN=1)
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import axi

PHI = (1 + 5 ** 0.5) / 2
R = json.load(open("out/axi_risultati.json"))
COL = {"A": "#2471a3", "B": "#17a589", "C": "#8e44ad", "D": "#d35400"}

fig = plt.figure(figsize=(13, 10), layout="constrained")
gs = fig.add_gridspec(3, 2, height_ratios=[1.15, 1, 1])

# 1) curve: migliore resistenza ottenibile vs L/D, normalizzata all'ottimo
ax = fig.add_subplot(gs[0, :])
for k, v in R.items():
    f = np.array([c["f"] for c in v["curva"]]); c = np.array([c["c"] for c in v["curva"]])
    rel = 100 * (c / v["ottimo"]["c"] - 1)
    ax.plot(f, rel, "o-", c=COL[k[0]], label=v["descr"], ms=4)
    ax.plot(v["ottimo"]["f"], 0, "*", c=COL[k[0]], ms=14, mec="k")
ax.axvline(PHI, ls="--", c="#c0392b"); ax.text(PHI + 0.05, ax.get_ylim()[1] * 0.92, "φ = 1,618", color="#c0392b")
ax.set_xscale("log"); ax.set_xticks([1.3, PHI, 2, 3, 4, 5, 6, 8, 10, 12])
ax.set_xticklabels(["1,3", "1,618", "2", "3", "4", "5", "6", "8", "10", "12"])
ax.set_ylim(-3, 120)
ax.set_xlabel("snellezza L/D"); ax.set_ylabel("resistenza in più rispetto all'ottimo [%]")
ax.set_title("Migliore resistenza ottenibile a ogni snellezza (★ = ottimo globale)", loc="left", fontsize=10)
ax.legend(fontsize=8, loc="upper right"); ax.grid(alpha=.3)
ax.spines[["top", "right"]].set_visible(False)

# 2) profili ottimi (volume fisso: stessa scala)
ax = fig.add_subplot(gs[1, :])
for k, v in R.items():
    x = np.array(v["ottimo"]["x"]); r = np.array(v["ottimo"]["r"])
    if v["obiettivo"] == "frontale":           # riscala a diametro uguale alla goccia A
        continue
    ax.plot(x, r, c=COL[k[0]], lw=2, label=f"{k[0]}: L/D {v['ottimo']['f']:.1f}, spessore max a {100*v['ottimo']['x_spess_max']:.0f}%")
    ax.plot(x, -r, c=COL[k[0]], lw=2)
xg, rg = axi.ellipsoid(PHI); ax.plot(xg, rg, c="#c0392b", ls="--", label="ellissoide aureo (stesso volume)"); ax.plot(xg, -rg, c="#c0392b", ls="--")
ax.set_aspect("equal"); ax.legend(fontsize=8, loc="upper right")
ax.set_title("Forme ottime a PARI VOLUME (sezione del solido di rivoluzione, flusso da sinistra)", loc="left", fontsize=10)
ax.spines[["top", "right"]].set_visible(False)

# 3) profilo ottimo a sezione frontale fissa
ax = fig.add_subplot(gs[2, 0])
v = R.get("D_frontale_ReD1e6")
if v:
    x = np.array(v["ottimo"]["x"]); r = np.array(v["ottimo"]["r"]); s = 0.5 / r.max()
    ax.plot(x * s, r * s, c=COL["D"], lw=2, label=f"ottimo D: L/D {v['ottimo']['f']:.2f}")
    ax.plot(x * s, -r * s, c=COL["D"], lw=2)
    xg, rg = axi.ellipsoid(PHI); s = 0.5 / rg.max()
    ax.plot(xg * s, rg * s, c="#c0392b", ls="--", label="ellissoide aureo"); ax.plot(xg * s, -rg * s, c="#c0392b", ls="--")
    ax.set_aspect("equal"); ax.legend(fontsize=8)
    ax.set_title("Pari SEZIONE FRONTALE (es. veicolo di altezza data)", loc="left", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)

# 4) tabella confronto
ax = fig.add_subplot(gs[2, 1]); ax.axis("off")
rows = []
for k, v in R.items():
    o = v["ottimo"]["c"]
    cf = v["confronto"]
    rows.append([k[0], f"{o:.4f}",
                 f"+{100*(cf['goccia aurea (L/D 1,618)']['c']/o-1):.0f}%" + ("*" if cf['goccia aurea (L/D 1,618)']['staccato'] else ""),
                 f"+{100*(cf['ellissoide aureo (L/D 1,618)']['c']/o-1):.0f}%" + ("*" if cf['ellissoide aureo (L/D 1,618)']['staccato'] else ""),
                 f"+{100*(cf['goccia standard L/D 4,5']['c']/o-1):.0f}%"])
t = ax.table(cellText=rows, colLabels=["scen.", "ottimo", "goccia\naurea", "ellissoide\naureo", "goccia\nL/D 4,5"],
             loc="center", cellLoc="center")
t.auto_set_font_size(False); t.set_fontsize(9); t.scale(1, 1.6)
ax.set_title("Resistenza in più rispetto all'ottimo\n(* flusso staccato: valore stimato)", fontsize=10)

fig.suptitle("Corpi di rivoluzione 3D a Reynolds reali — pannelli + strato limite integrale", fontsize=12)
# layout="constrained" above: no overlap between the axis label and the next panel title
plt.savefig("out/ottimo_3D.png", dpi=110)
print("ok")
