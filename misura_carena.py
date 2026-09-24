"""Perdita d'imbocco K e bilancio di potenza: propulsore intubato vs disco libero.

Per ciascun caso (campi in <caso>/risultati/{C,U,p}):
  sezione del ventilatore (appena a monte del disco): portata Q e pressione totale
  media pesata sulla portata p0_vent;
  stazione a monte della presa: il tubo di flusso più vicino alla parete che porta
  la stessa portata Q (quello che entra nel condotto): p0_cattura;
  K = (p0_cattura - p0_vent) / (rho U^2 / 2).
Potenza in autopropulsione: P = T * <u>_disco (per cuneo di 5 gradi); il confronto
è con il corpo senza propulsore spinto da un'elica ideale: P_conv = D_nudo * U.

uso: python misura_carena.py
"""
import json, sys
import numpy as np
from scipy.interpolate import griddata
from plot_cfd import campo

P0INF = 0.5
D_NUDO = 6.9574e-06                    # corpo nudo (conv_ref, stessa turbolenza di fondo), per cuneo
X_CATTURA = 0.835


def carica(caso):
    C = campo(f"{caso}/risultati/C", True); U = campo(f"{caso}/risultati/U", True)
    p = campo(f"{caso}/risultati/p")
    return C[:, 0], np.hypot(C[:, 1], C[:, 2]), U, p


def sezione(dati, xs, r0, r1, n=600, fin=0.003):
    """u_x e p0 lungo la linea x = xs tra i raggi r0 (parete) e r1; la finestra in x
    si allarga finché contiene celle da entrambi i lati."""
    x, r, U, p = dati
    while True:
        m = (np.abs(x - xs) < fin) & (r > r0 - 0.003) & (r < r1 + 0.003)
        if (x[m] < xs).any() and (x[m] > xs).any():
            break
        fin *= 1.5
    sc = 0.02                                                # comprime x: triangoli meno schiacciati
    pts = np.c_[(x[m] - xs) * sc, r[m]]
    rr = r0 + (r1 - r0) * (1 - np.cos(np.linspace(0, np.pi, n))) / 2
    q = np.c_[np.zeros(n), rr]
    ux = griddata(pts, U[m, 0], q, method="linear")
    p0 = griddata(pts, p[m] + 0.5 * (U[m] ** 2).sum(1), q, method="linear")
    p0n = griddata(pts, p[m] + 0.5 * (U[m] ** 2).sum(1), q, method="nearest")
    ux = np.where(np.isnan(ux), 0.0, ux); p0 = np.where(np.isnan(p0), p0n, p0)
    return rr, ux, p0


def integrali(rr, ux, p0):
    dQ = 2 * np.pi * rr * ux
    Q = np.trapezoid(dQ, rr)
    return Q, np.trapezoid(dQ * p0, rr) / Q


def cattura(dati, rb, Q):
    rr, ux, p0 = sezione(dati, X_CATTURA, rb, rb + 0.25, n=3000)
    dQ = 2 * np.pi * rr * ux
    Qc = np.concatenate([[0], np.cumsum(0.5 * (dQ[1:] + dQ[:-1]) * np.diff(rr))])
    k = int(np.searchsorted(Qc, Q))
    r_cap = float(np.interp(Q, Qc[k - 1:k + 1], rr[k - 1:k + 1]))
    s = rr <= r_cap
    Qx, p0c = integrali(np.r_[rr[s], r_cap], np.r_[ux[s], np.interp(r_cap, rr, ux)],
                        np.r_[p0[s], np.interp(r_cap, rr, p0)])
    return r_cap, p0c


def autopropulsione(fn):
    righe = [l for l in open(fn) if l.startswith("iter")]
    ult = righe[-1]
    T = float(ult.split("T=")[1].split()[0]); D = float(ult.split("D=")[1].split()[0])
    u = float(ult.split("=")[-1])
    return T, D, u


def analizza(nome, caso, xf, r_esterno, rb_fun):
    dati = carica(caso)
    rr, ux, p0 = sezione(dati, xf, rb_fun(xf), r_esterno)
    Q, p0v = integrali(rr, ux, p0)
    r_cap, p0c = cattura(dati, rb_fun(X_CATTURA), Q)
    T, D, u = autopropulsione(f"{caso}/autopropulsione.txt")
    P = T * u
    out = dict(caso=nome, Q=Q, p0_cattura=p0c, p0_ventilatore=p0v,
               K=(p0c - p0v) / 0.5, r_cattura=r_cap, spessore_catturato=r_cap - rb_fun(X_CATTURA),
               T=T, D=D, u_disco=u, P=P, P_conv=D_NUDO, risparmio=100 * (1 - P / D_NUDO),
               velocita_media_vent=Q / np.trapezoid(2 * np.pi * rr, rr))
    print(f"{nome:28s} Q={out['Q']:.3e}  p0 cattura {p0c:+.4f} -> ventilatore {p0v:+.4f}  "
          f"K = {100 * out['K']:5.1f}%  | T={T:.4e} <u>={u:.4f}  P/P_conv = {P / D_NUDO:.4f} "
          f"({out['risparmio']:+.1f}% di potenza risparmiata)")
    return out


if __name__ == "__main__":
    prof = json.load(open("out/cfd_profili.json"))["conv_1e6f"]
    xb = np.array(prof["x"]) / prof["x"][-1]; rbb = np.array(prof["r"]) / prof["x"][-1]
    rb = lambda xx: float(np.interp(xx, xb, rbb))
    res = []
    # disco libero (senza carenatura): disco in x 0.895-0.905, raggio 0.03
    res.append(analizza("disco libero (R=0.03)", "cfd/disco_A", 0.891, 0.03, rb))
    caso = sys.argv[1] if len(sys.argv) > 1 else "cfd/carena_D"
    g = json.load(open(f"{caso}/carena.json"))
    import cfd_carena as cc
    for k in ("X_LE", "C", "TAU", "TTE", "R_LE", "R_TE"):       # geometria del caso misurato
        setattr(cc, k, g[k])
    xf = g["X_F"] - g["SP_F"] / 2 - 0.004
    res.append(analizza("presa intubata con labbro", caso, xf, float(cc.carena(xf)[0]), rb))
    json.dump(res, open("out/misura_carena.json", "w"), indent=1)
