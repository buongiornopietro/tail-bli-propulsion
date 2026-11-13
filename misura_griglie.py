"""Nota 1: carena D sulle tre griglie (rada, base, fine), ognuna confrontata con il SUO corpo nudo.
Stessa misura di misura_v3.py (misura_carena.analizza) e stesse grandezze navali di rivista_navale.py
(identita' di spinta), ma su cartelle date: non tocca out/misura_v3.json ne' out/rivista_navale.json.

uso: python misura_griglie.py nome cartella D_nudo [nome cartella D_nudo ...]   -> out/nota1_griglie.json
"""
import json, os, sys
import numpy as np
import misura_carena as mc
import cfd_carena as cc
from confronto_elica import fattore_disco

SPICCHI = 72                       # come rivista_navale.py: forze per spicchio di 5 gradi -> anello intero
prof = json.load(open("out/cfd_profili.json"))["conv_1e6f"]
xb = np.array(prof["x"]) / prof["x"][-1]; rbb = np.array(prof["r"]) / prof["x"][-1]
rb = lambda xx: float(np.interp(xx, xb, rbb))


def carena(nome, caso, D_nudo):    # come misura_v3.carena
    g = json.load(open(f"{caso}/carena.json"))
    for k in ("X_LE", "C", "TAU", "TTE", "R_LE", "R_TE"):
        setattr(cc, k, g[k])
    mc.D_NUDO = D_nudo
    xf = g["X_F"] - g["SP_F"] / 2 - 0.004
    return mc.analizza(nome, caso, xf, float(cc.carena(xf)[0]), rb)


a = sys.argv[1:]
uscita = "out/nota1_griglie.json"
ris = json.load(open(uscita)) if os.path.exists(uscita) else {}
for nome, caso, D0 in zip(a[0::3], a[1::3], map(float, a[2::3])):
    r = carena(nome, caso, D0)
    A = r["Q"] / r["velocita_media_vent"]
    f = fattore_disco(r["P_conv"], A)[1]
    T, u = r["T"], r["u_disco"]
    t, w_tot = 1 - D0 / T, 1 - u
    V_A = u - SPICCHI * T / (2 * A * u)
    w_T = 1 - V_A
    eta_H, eta_O = (1 - t) / (1 - w_T), V_A / u
    r.update(cartella=caso, A=A, variazione_ideale=100 * (r["P"] / r["P_conv"] - 1),
             variazione_pari_area=100 * (r["P"] / (r["P_conv"] * f) - 1),
             t=t, w_tot=w_tot, w_T=w_T, eta_H=eta_H, eta_O=eta_O, eta_D=eta_H * eta_O)
    ris[nome] = r
    print(f"{nome:10s} T {T:.6e} D {r['D']:.6e} u {u:.4f} K {100 * r['K']:5.2f}%  P vs ideale "
          f"{r['variazione_ideale']:+6.2f}%  vs pari area {r['variazione_pari_area']:+6.2f}%  t {t:+.4f}  "
          f"w_T {w_T:+.4f}  eta_H {eta_H:.4f}  eta_O {eta_O:.4f}  eta_D {eta_H * eta_O:.4f}")
json.dump(ris, open(uscita, "w"), indent=1)
