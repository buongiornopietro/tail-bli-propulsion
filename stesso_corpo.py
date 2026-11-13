"""Revisione Fable xhigh, 12/11/2026 (F5); adattato ed eseguito dalla Regia col via di Pietro (12/11 10:05 UTC).

Soglia K* "a parita' di corpo": propulsore in coda applicato al corpo convenzionale ottimo (forma fissa),
ottimizzando solo la posizione della presa x_s e la portata Q. Da confrontare con l'ottimo congiunto dei JSON v4
(forma libera), che finisce con una coda tronca non verificata dal modello (PER_LA_REGIA.md, sezione F5).

A forma, x_s e Q fissi C_P e' lineare in K (D0): basta la mappa CP0(x_s, Q) a K = 0, poi
CP(K) = min su (x_s, Q) di CP0 + K Q / eta.

Costo stimato: 12 casi (6 Re x 2 ipotesi d'ingestione), circa 1000 valutazioni del modello ciascuno
(mappa 30 x 30 piu' rifinitura), circa 0,1 s l'una: una decina di minuti con 4 processi. Stima, non misura.

Versione della Regia: bli.py ufficiale (correzioni di E e di Young sempre attive) e PHI_MODE.
(La copia della revisione con FIX_E, FIX_D riproduce i C_P dei JSON v4 a sei cifre, riproduci.py.)

Uso: python stesso_corpo.py <cartella dei JSON v4> <cartella con axi.py, bli.py, bli_opt.py> <uscita.json>
"""
import json, sys, os, time
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from multiprocessing import Pool

ETA = 0.85
KS = [0.0, 0.025, 0.05, 0.075, 0.10, 0.15, 0.20]
RES = [3e5, 1e6, 3e6, 1e7, 3e7, 1e8]


def lavoro(a):
    Re, modo, D, M = a
    sys.path.insert(0, M)
    import bli, bli_opt
    from scipy.optimize import minimize
    bli.PHI_MODE = modo                       # bli.py ufficiale: correzioni di E e di Young sempre attive
    d = json.load(open(os.path.join(D, f'soglia_Re{Re:.0e}_forzata_v4{"lin" if modo == "lineare" else "pro"}.json')))
    n1 = 1 + bli_opt.NW
    pc = np.clip(d['p_conv'], bli_opt.LO[:n1], bli_opt.HI[:n1]); cc = d['C_conv']

    def cp0(xs, lq, det=False):
        x, r, xs_, Q = bli_opt.unpack(np.concatenate([pc, [xs, lq]]))
        try:
            return bli.power(x, r, Re, 'forzata', xs_, Q, detail=det, eta_pompa=ETA, K_condotto=0.0)
        except Exception:
            return None if det else 10.0

    XS = np.linspace(0.55, 0.96, 30); LQ = np.linspace(np.log(0.01), np.log(0.25), 30)
    G = np.array([[cp0(a_, b_) for b_ in LQ] for a_ in XS])
    ris = []
    for K in KS:
        T = G + K * np.exp(LQ)[None, :] / ETA
        i, j = np.unravel_index(np.argmin(T), T.shape)
        cl = lambda q: (float(np.clip(q[0], .55, .96)), float(np.clip(q[1], LQ[0], LQ[-1])))
        f = lambda q: cp0(*cl(q)) + K * np.exp(cl(q)[1]) / ETA
        l = minimize(f, [XS[i], LQ[j]], method='Nelder-Mead', options=dict(maxfev=120, xatol=1e-4, fatol=1e-8))
        xs, lq = cl(l.x) if l.fun < T[i, j] else (XS[i], LQ[j])
        dd = cp0(xs, lq, det=True)
        c = dd['CP'] + K * np.exp(lq) / ETA
        ris.append(dict(K=K, CP=float(c), guadagno=float(100 * (c / cc - 1)), xs=float(xs), Q=float(np.exp(lq)),
                        Q_BL=dd['Q_BL'], phi=dd['phi'], pen=dd['penalita'], Vj=dd['Vj']))
    # K* esatto (nota della revisione): a corpo fisso CP(K) = min su (x_s, Q) di [CP0 + K Q/eta], concava in K; zero
    # sulla mappa, poi rifinitura del progetto allo zero: K* = zero della sua retta (l'inviluppo puo' solo salire).
    from scipy.optimize import brentq
    U = lambda k: float(np.min(G + k * np.exp(LQ)[None, :] / ETA)) - cc
    k_map = brentq(U, 0, 0.6) if U(0) < 0 < U(0.6) else float('nan')
    k_star = k_map
    if np.isfinite(k_map):
        T = G + k_map * np.exp(LQ)[None, :] / ETA
        i, j = np.unravel_index(np.argmin(T), T.shape)
        cl = lambda q: (float(np.clip(q[0], .55, .96)), float(np.clip(q[1], LQ[0], LQ[-1])))
        f = lambda q: cp0(*cl(q)) + k_map * np.exp(cl(q)[1]) / ETA
        l = minimize(f, [XS[i], LQ[j]], method='Nelder-Mead', options=dict(maxfev=200, xatol=1e-5, fatol=1e-9))
        xs, lq = cl(l.x)
        c0 = cp0(xs, lq)
        k_star = max(k_map, (cc - c0) / (np.exp(lq) / ETA))
        dd = cp0(xs, lq, det=True)
        stella = dict(K_star=float(k_star), K_mappa=float(k_map), xs=float(xs), Q=float(np.exp(lq)), Q_BL=dd['Q_BL'],
                      phi=dd['phi'], pen=dd['penalita'])
    else:
        stella = dict(K_star=k_star)
    return dict(Re=Re, modo=modo, C_conv=cc, stesso_corpo=ris, stella=stella,
                congiunto=[dict(K=e['K'], CP=e['CP'], Q=e['Q'], xs=e['xs']) for e in d['bli']],
                mappa=dict(XS=XS.tolist(), LQ=LQ.tolist(), G=G.tolist()))


if __name__ == '__main__':
    D, M, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
    t0 = time.time()
    with Pool(4) as po:
        R = po.map(lavoro, [(Re, m, D, M) for m in ('lineare', 'profilo') for Re in RES], chunksize=1)
    json.dump(R, open(OUT, 'w'), indent=1)
    print('fatto in', round(time.time() - t0), 's')
