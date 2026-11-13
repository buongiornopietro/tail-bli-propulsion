"""Figura 4 della nota 1, versione v4 (12/11/2026, revisione Fable xhigh sezione H3).
Fascia piena: K* sullo stesso corpo convenzionale ottimo (stesso_corpo.py -> out/stesso_corpo_v4.json), bordo basso
"whole boundary layer ingested", bordo alto "inner layer only". Fascia tratteggiata senza riempimento: corpo
riottimizzato insieme al propulsore (out/fascia_K_v4.json, regola della revisione), stima superiore condizionata;
a Re_V 1e8, strato intero, solo la famiglia snella (semi 0 e 2: 3,3). Etichette in inglese (figura della rivista).
Scrive out/mappa_finale_v4_en.png."""
import json
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

RES = [3e5, 1e6, 3e6, 1e7, 3e7, 1e8]
S = {(r["Re"], r["modo"]): 100 * r["stella"]["K_star"] for r in json.load(open("out/stesso_corpo_v4.json"))}
F = json.load(open("out/fascia_K_v4.json"))
s_whole = [S[(r, "lineare")] for r in RES]; s_inner = [S[(r, "profilo")] for r in RES]
j_whole = [F[f"lin {r:.0e}"]["valore"] for r in RES]; j_inner = [F[f"pro {r:.0e}"]["valore"] for r in RES]
j_whole[-1] = 3.3          # 1e8: famiglia snella (semi 0 e 2), dichiarato in didascalia

fig, ax = plt.subplots(figsize=(8.5, 5.0))
ax.fill_between(RES, s_whole, s_inner, color="#d35400", alpha=.25, lw=0)
ax.plot(RES, s_whole, "o-", c="#d35400", lw=2.2, ms=6, label="same conventional body: whole boundary layer ingested")
ax.plot(RES, s_inner, "s-", c="#a04000", lw=1.6, ms=5, label="same conventional body: inner layer only")
ax.plot(RES, j_whole, "o--", c="#5d6d7e", lw=1.3, ms=5, mfc="none",
        label="re-optimised body: whole boundary layer ingested")
ax.plot(RES, j_inner, "s--", c="#5d6d7e", lw=1.3, ms=5, mfc="none",
        label="re-optimised body: inner layer only")
ax.text(0.99, 0.76, "re-optimised body: truncated afterbody,\nnot verified by the model (upper estimate)",
        transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color="#5d6d7e")
ax.set_xscale("log"); ax.set_xlim(2e5, 1.5e8); ax.set_ylim(0, 13)
ax.set_xlabel("volumetric Reynolds number  Re_V")
ax.set_ylabel("break-even duct loss  K*  [% of ½ρU²]")
ax.set_title("Tail BLI pays below the curves, the conventional propeller above them", fontsize=10, loc="left")
ax.legend(fontsize=7.5, loc="upper right"); ax.grid(alpha=.3)
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout(); plt.savefig("out/mappa_finale_v4_en.png", dpi=150)
for n, v in (("same body, whole", s_whole), ("same body, inner", s_inner), ("re-optimised, whole", j_whole),
             ("re-optimised, inner", j_inner)):
    print(n, [f"{x:.1f}" for x in v])
