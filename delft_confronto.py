"""Confronto con l'esperimento di Delft (Della Corte et al., AIAA J. 59(12), 2021).

Grandezza confrontabile: il flusso di energia cinetica depositata nella scia (termini
Ė_a + Ė_v del bilancio di potenza di Drela 2009), riferito alla potenza D*U, sul piano
del bordo d'uscita e poco a valle. Delft lo chiama C_E e misura C_E / C_Fx circa 0.13 sul
corpo nudo (L/D circa 11, PIV a 20 m/s cioè Re_L circa 2.4e6, campagna 1.7-4.7e6, innesco al 5%):
"circa il 13% delle perdite totali di potenza avviene nella scia", cioè il beneficio
ideale del BLI per il solo corpo.

Qui, dal RANS del corpo nudo (campi esportati in cfd/<caso>/risultati, nessuna
simulazione):
  Ė_a = ∫ ½ u (u - U)² dA ,  Ė_v = ∫ ½ u (v_r² + v_t²) dA ,  su r da parete (o asse) a R_MAX
  Ė_p = ∫ (p - p_inf)(u - U) dA (lavoro di pressione): il C_E di Delft (eq. 10) è Ė_a+Ė_v+Ė_p,
  riferito al flusso di deficit di quantità di moto F_x = ∫ u (U - u) dA - ∫ p dA alla stessa stazione.

uso: python delft_confronto.py [caso] [D_per_cuneo]
"""
import json, sys
import numpy as np
from scipy.interpolate import griddata
import misura_carena as mc

caso = sys.argv[1] if len(sys.argv) > 1 else "cfd/conv_ref"
D_CUNEO = float(sys.argv[2]) if len(sys.argv) > 2 else 6.9574e-06
CUNEI = 360 / 5
R_MAX = 0.25

prof = json.load(open("out/cfd_profili.json"))["conv_1e6f"]
xb = np.array(prof["x"]) / prof["x"][-1]; rbb = np.array(prof["r"]) / prof["x"][-1]
rb = lambda xx: float(np.interp(xx, xb, rbb)) if xx <= 1 else 0.0

C = mc.campo(f"{caso}/risultati/C", True); U = mc.campo(f"{caso}/risultati/U", True)
p = mc.campo(f"{caso}/risultati/p")
x, r = C[:, 0], np.hypot(C[:, 1], C[:, 2])
ur = (C[:, 1] * U[:, 1] + C[:, 2] * U[:, 2]) / np.maximum(r, 1e-12)
ut = (C[:, 1] * U[:, 2] - C[:, 2] * U[:, 1]) / np.maximum(r, 1e-12)


def piano(xs, r0, n=2000, fin=0.003):
    while True:
        m = (np.abs(x - xs) < fin) & (r < R_MAX + 0.01)
        if (x[m] < xs).any() and (x[m] > xs).any():
            break
        fin *= 1.5
    pts = np.c_[(x[m] - xs) * 0.02, r[m]]
    rr = r0 + (R_MAX - r0) * (1 - np.cos(np.linspace(0, np.pi, n))) / 2
    q = np.c_[np.zeros(n), rr]
    f = lambda v: np.nan_to_num(griddata(pts, v[m], q, method="linear"),
                                nan=0.0)
    return rr, f(U[:, 0]), f(ur), f(ut), f(p)


D = CUNEI * D_CUNEO
out = []
for xs in (0.95, 1.0, 1.02, 1.05, 1.1, 1.2, 1.4):
    rr, u, vr, vt, pp = piano(xs, rb(xs))
    dA = 2 * np.pi * rr
    Ea = np.trapezoid(dA * 0.5 * u * (u - 1) ** 2, rr)
    Ev = np.trapezoid(dA * 0.5 * u * (vr ** 2 + vt ** 2), rr)
    Ep = np.trapezoid(dA * pp * (u - 1), rr)          # lavoro di pressione (Drela 2009), p_inf = 0
    mom = np.trapezoid(dA * (u * (1 - u) - pp), rr)
    E = Ea + Ev + Ep                                   # come C_E di Delft (eq. 10: a + r + t + p)
    out.append(dict(x=xs, Ea=float(Ea), Ev=float(Ev), Ep=float(Ep), frazione=float((Ea + Ev) / D),
                    frazione_con_Ep=float(E / D), frazione_su_qdm=float(E / mom),
                    quantita_di_moto_su_D=float(mom / D)))
    print(f"x = {xs:4.2f}   (Ea+Ev)/(D U) = {100 * (Ea + Ev) / D:5.1f}%   Ep/(D U) = {100 * Ep / D:5.1f}%   "
          f"(Ea+Ev+Ep)/(D U) = {100 * E / D:5.1f}%   /(F_x U) = {100 * E / mom:5.1f}%   q.d.m./D = {mom / D:5.2f}")
json.dump(dict(caso=caso, D=D, stazioni=out, delft_CE_su_CFx=0.13),
          open("out/delft_confronto.json", "w"), indent=1)
