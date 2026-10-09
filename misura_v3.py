"""Versione 3 della nota: K e potenza in autopropulsione dei casi v3 (calcolati sull'EVO-X2),
con le stesse funzioni di misura_carena.py; potenza riferita al corpo nudo calcolato con la
stessa griglia e lo stesso Reynolds, sia con l'elica ideale (P = D U) sia con un'elica isolata
della stessa area del propulsore (teoria del disco, come confronto_elica.py).
Scrive out/misura_v3.json. Nessuna simulazione.

uso: python misura_v3.py
"""
import json
import numpy as np
import misura_carena as mc
import cfd_carena as cc
from confronto_elica import fattore_disco

D_BASE = 6.9574e-06          # conv_ref, griglia base
D_RADA = 6.93443e-06         # v3_conv_rada, 2400 it. (variazione < 1e-4 relativa nelle ultime 1000)
D_RE7 = 4.9111e-06           # v3_conv_re7, 2400 it.

prof = json.load(open("out/cfd_profili.json"))["conv_1e6f"]
xb = np.array(prof["x"]) / prof["x"][-1]; rbb = np.array(prof["r"]) / prof["x"][-1]
rb = lambda xx: float(np.interp(xx, xb, rbb))


def carena(nome, caso, D_nudo):
    g = json.load(open(f"{caso}/carena.json"))
    for k in ("X_LE", "C", "TAU", "TTE", "R_LE", "R_TE"):
        setattr(cc, k, g[k])
    mc.D_NUDO = D_nudo
    xf = g["X_F"] - g["SP_F"] / 2 - 0.004
    return mc.analizza(nome, caso, xf, float(cc.carena(xf)[0]), rb)


def disco(nome, caso, D_nudo):
    mc.D_NUDO = D_nudo
    return mc.analizza(nome, caso, 0.891, 0.03, rb)


ris = [carena("carenatura E (corta e sottile)", "cfd/v3_carena_E", D_BASE),
       carena("carenatura F (meno convergente)", "cfd/v3_carena_F", D_BASE),
       carena("carenatura G (quasi a sezione costante)", "cfd/v3_carena_G", D_BASE),
       disco("disco libero, griglia rada", "cfd/v3_disco_rada", D_RADA),
       disco("disco libero, Re x10", "cfd/v3_disco_re7", D_RE7),
       carena("carenatura D, Re x10 (NON convergente)", "cfd/v3_carena_re7", D_RE7)]
print()
for r in ris:
    A = r["Q"] / r["velocita_media_vent"]
    f = fattore_disco(r["P_conv"], A)[1]
    r["A"] = A
    r["variazione_ideale"] = 100 * (r["P"] / r["P_conv"] - 1)
    r["variazione_pari_area"] = 100 * (r["P"] / (r["P_conv"] * f) - 1)
    print(f"{r['caso']:42s} K {100 * r['K']:5.1f}%  u {r['u_disco']:.3f}  T/D_nudo {r['T'] / r['P_conv']:.3f}  "
          f"P vs ideale {r['variazione_ideale']:+5.1f}%  vs pari area {r['variazione_pari_area']:+5.1f}%")
json.dump(ris, open("out/misura_v3.json", "w"), indent=1)
