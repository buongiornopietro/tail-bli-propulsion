"""Verifica di griglia: K e potenza in autopropulsione sui casi fini, riferiti al corpo nudo
calcolato sulla stessa griglia (conv_fine). Usa le funzioni di misura_carena.py.
uso: python misura_fine.py [caso_carena] [caso_disco] [D_nudo]"""
import json, sys
import numpy as np
import misura_carena as mc

arg = sys.argv[1:] + [None] * 3
caso = arg[0] or "cfd/carena_fine2"
mc.D_NUDO = float(arg[2] or 6.97255e-06)               # conv_fine, 2400 iterazioni
prof = json.load(open("out/cfd_profili.json"))["conv_1e6f"]
xb = np.array(prof["x"]) / prof["x"][-1]; rbb = np.array(prof["r"]) / prof["x"][-1]
rb = lambda xx: float(np.interp(xx, xb, rbb))
res = [mc.analizza("disco libero (fine)", arg[1] or "cfd/disco_fine", 0.891, 0.03, rb)]
g = json.load(open(f"{caso}/carena.json"))
import cfd_carena as cc
for k in ("X_LE", "C", "TAU", "TTE", "R_LE", "R_TE"):
    setattr(cc, k, g[k])
xf = g["X_F"] - g["SP_F"] / 2 - 0.004
res.append(mc.analizza("presa intubata (fine)", caso, xf, float(cc.carena(xf)[0]), rb))
json.dump(res, open("out/misura_fine.json", "w"), indent=1)
