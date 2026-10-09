"""y+ della prima cella sul corpo, stimato dai campi esportati (solo lettura, nessun calcolo).
Per ogni striscia in x prende la cella più vicina alla parete: distanza y = r_c - r_corpo(x_c)
(approssimazione normale-radiale), velocità |U|; nel sottostrato viscoso u+ = y+, quindi
y+ = sqrt(|U| * y * Re_L).
uso: python yplus_campi.py caso Re_L
"""
import json, sys
import numpy as np
import misura_carena as mc

caso, Re = sys.argv[1], float(sys.argv[2])
prof = json.load(open("out/cfd_profili.json"))["conv_1e6f"]
xb = np.array(prof["x"]) / prof["x"][-1]; rbb = np.array(prof["r"]) / prof["x"][-1]
x, r, U, p = mc.carica(caso)
m = (x > 0.02) & (x < 0.98)
y = r[m] - np.interp(x[m], xb, rbb); um = np.linalg.norm(U[m], axis=1); xm = x[m]
ok = y > 0
y, um, xm = y[ok], um[ok], xm[ok]
bins = np.linspace(0.02, 0.98, 97)
yp = []
for a, b in zip(bins[:-1], bins[1:]):
    s = (xm >= a) & (xm < b)
    if s.any():
        i = np.argmin(y[s])
        yp.append(np.sqrt(um[s][i] * y[s][i] * Re))
yp = np.array(yp)
print(f"{caso}: y+ prima cella (0.02<x<0.98)  media {yp.mean():.2f}  max {yp.max():.2f}")
