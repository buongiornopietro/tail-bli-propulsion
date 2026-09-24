"""Misura della perdita d'imbocco K di una presa che ingerisce lo strato limite.

K = (Delta_p0_misurata - Delta_p0_ideale) / (rho U^2 / 2)

  Delta_p0_misurata = P_pompa / Q          (dalla CFD: flussi di pressione totale)
  Delta_p0_ideale   = (1/Q) int (p0_inf - p0) u dA   sulla porzione di strato
                      limite che porta la stessa portata Q, misurata a monte
                      della presa nello stesso campo CFD.

uso: python misura_K.py <caso1> [<caso2> ...]
     ogni caso deve avere risultati/{C,U,p} e bilanci.json
"""
import json, sys
import numpy as np
from plot_cfd import campo

CUNEO = 5.0 / 360.0
P0INF = 0.5                      # p_inf = 0, U = 1 (grandezze cinematiche)


def profilo(d, x_staz, x_corpo, r_corpo, dx=0.002):
    """Profilo dello strato limite: (r, u_x, p0) ordinati dalla parete."""
    C = campo(f"{d}/C", True); U = campo(f"{d}/U", True); p = campo(f"{d}/p")
    x = C[:, 0]; r = np.hypot(C[:, 1], C[:, 2])
    rb = float(np.interp(x_staz, x_corpo, r_corpo))
    # distanza dalla parete calcolata cella per cella (il raggio del corpo varia
    # lungo x: usare il raggio della sola stazione mescolerebbe punti diversi)
    m = np.abs(x - x_staz) < dx
    dist = r[m] - np.interp(x[m], x_corpo, r_corpo)
    ok = (dist > 0) & (dist < 0.25)
    o = np.argsort(dist[ok])
    return rb, rb + dist[ok][o], U[m][ok][o][:, 0], p[m][ok][o] + 0.5 * (U[m][ok][o] ** 2).sum(1)


def misura(caso):
    info = json.load(open(f"{caso}/gold.json"))
    bil = json.load(open(f"{caso}/bilanci.json"))
    prof = json.load(open("out/cfd_profili.json"))[bil["profilo"]]
    xb = np.array(prof["x"]); rb_ = np.array(prof["r"]); L = xb[-1]
    x_corpo, r_corpo = xb / L, rb_ / L
    Q = info["Q_L"]                                   # portata (corpo intero)
    P_cuneo = -(bil["p0_presa"] + bil["p0_scarico"])  # potenza pompa per cuneo
    P = P_cuneo / CUNEO                               # corpo intero
    dp0_mis = P / Q

    x_staz = info["xa"] - 0.05
    rb, rr, ux, p0 = profilo(caso + "/risultati", x_staz, x_corpo, r_corpo)
    dA = 2 * np.pi * rr * np.gradient(rr)
    dQ = ux * dA
    Qcum = np.cumsum(dQ)
    if Qcum[-1] < Q:
        raise SystemExit(f"{caso}: strato limite troppo sottile per la portata richiesta")
    k = int(np.searchsorted(Qcum, Q))
    w = np.ones(k + 1); w[-1] = (Q - (Qcum[k - 1] if k else 0)) / dQ[k]
    dp0_id = float(np.sum((P0INF - p0[:k + 1]) * dQ[:k + 1] * w) / Q)
    # pressione totale media (pesata sulla portata) sulla presa
    p0_presa = bil["p0_presa"] / bil["Q_presa"]
    p0_monte = P0INF - dp0_id
    K = (p0_monte - p0_presa) / 0.5              # perdita d'imbocco
    getto = (bil["p0_scarico"] / bil["Q_scarico"] - P0INF) / 0.5   # eccesso del getto
    return dict(caso=caso, larghezza=info["xb"] - info["xa"], p0_monte=p0_monte,
                p0_presa=p0_presa, eccesso_getto=getto,
                v_ingresso=Q / (2 * np.pi * np.interp(info["xa"], x_corpo, r_corpo) * (info["xb"] - info["xa"])),
                Q=Q, dp0_misurata=dp0_mis, dp0_ideale=dp0_id, K=K,
                spessore_ingerito=float(rr[k] - rb))


if __name__ == "__main__":
    out = []
    for c in sys.argv[1:]:
        r = misura(c)
        out.append(r)
        print(f"{c}: larghezza {r['larghezza']:.3f} L, ingresso {r['v_ingresso']:.3f} U | "
              f"p0 a monte {r['p0_monte']:+.3f} -> sulla presa {r['p0_presa']:+.3f}  =>  "
              f"K = {100 * r['K']:.1f}%   (getto: {100 * r['eccesso_getto']:+.0f}%, "
              f"strato ingerito {r['spessore_ingerito']:.4f} L)")
    json.dump(out, open("out/misura_K.json", "w"), indent=1)
