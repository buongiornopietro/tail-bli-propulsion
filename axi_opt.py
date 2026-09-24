"""Ottimizzazione di corpi di rivoluzione a Reynolds reali (vedi axi.py).

Parametri: ln(L/D) e ln dei 4 pesi CST (w_0 = 1) -> 5 variabili.
Per ogni scenario:
  - ottimo globale con evoluzione differenziale
  - curva "migliore resistenza ottenibile a L/D fissato" (Nelder-Mead sui pesi)
  - forme di confronto: ellissoide aureo, goccia aurea, goccia standard L/D 4,5
"""
import json, time, sys
import numpy as np
from scipy.optimize import differential_evolution, minimize
import axi

PHI = (1 + 5 ** 0.5) / 2
NW = 4
LB = np.log(0.1); UB = np.log(5.0)

SCENARI = {
    "A_turbolento_ReV1e7": dict(Re=1e7, transition="forzata", obiettivo="volume",
                                descr="Volume fisso, strato limite turbolento, Re_V = 10^7"),
    "B_naturale_ReV1e7":   dict(Re=1e7, transition="naturale", obiettivo="volume",
                                descr="Volume fisso, transizione naturale, Re_V = 10^7"),
    "C_naturale_ReV1e6":   dict(Re=1e6, transition="naturale", obiettivo="volume",
                                descr="Volume fisso, transizione naturale, Re_V = 10^6"),
    "D_frontale_ReD1e6":   dict(Re=1e6, transition="forzata", obiettivo="frontale",
                                descr="Sezione frontale fissa, turbolento, Re_D = 10^6"),
}


def cost_shape(x, r, sc):
    """Coefficiente da minimizzare per lo scenario (V = 1 in x, r)."""
    if sc["obiettivo"] == "volume":
        return axi.drag(x, r, sc["Re"], sc["transition"])
    D = 2 * r.max()
    CDV = axi.drag(x, r, sc["Re"] / D, sc["transition"])     # Re_V = Re_D V^(1/3)/D
    return axi.CD_frontal(CDV, x, r)


BOUNDS = [(np.log(1.3), np.log(12.0))] + [(LB, UB)] * NW


def set_nw(n):
    """Cambia il numero di pesi CST (flessibilità della forma)."""
    global NW, BOUNDS
    NW = n
    BOUNDS = [(np.log(1.3), np.log(12.0))] + [(LB, UB)] * NW


PENDENZA_MAX = None      # se impostata, penalizza |dr/dx| oltre il limite


def penalita_forma(x, r):
    if PENDENZA_MAX is None:
        return 0.0
    sl = np.abs(np.gradient(r, x))[5:-5]          # esclusi naso e coda
    ecc = np.maximum(sl - PENDENZA_MAX, 0).max()
    return 10.0 * ecc


def cost(p, sc):
    p = np.clip(p, [b[0] for b in BOUNDS], [b[1] for b in BOUNDS])[:1 + NW]   # limiti sempre
    try:
        x, r = axi.body([np.exp(p[0])] + list(np.exp(p[1:])))
        c = cost_shape(x, r, sc) + penalita_forma(x, r)
        return float(c) if np.isfinite(c) else 10.0
    except Exception:
        return 10.0


def run(name, sc):
    t0 = time.time()
    res = differential_evolution(cost, BOUNDS, args=(sc,), popsize=10, maxiter=40,
                                 tol=1e-4, seed=3, polish=False, init="sobol")
    p = res.x
    # rifinitura locale
    loc = minimize(cost, p, args=(sc,), method="Nelder-Mead",
                   options=dict(maxfev=300, xatol=1e-3, fatol=1e-6))
    if loc.fun < res.fun:
        p = loc.x
    p = np.clip(p, [b[0] for b in BOUNDS], [b[1] for b in BOUNDS])
    best = cost(p, sc)
    x, r = axi.body([np.exp(p[0])] + list(np.exp(p[1:])))
    Dre = sc["Re"] if sc["obiettivo"] == "volume" else sc["Re"] / (2 * r.max())
    det = axi.drag(x, r, Dre, sc["transition"], detail=True)
    i_max = int(np.argmax(r))
    print(f"[{name}] ottimo: coeff={best:.5f}  L/D={det['f']:.2f}  "
          f"spessore max a x/L={x[i_max] / x[-1]:.2f}  transizione x/L={det['x_trans']:.2f}  "
          f"distacco={det['x_sep_turb']}  ({time.time() - t0:.0f}s, {res.nfev} valutazioni)",
          flush=True)

    # curva: migliore resistenza a L/D fissato
    curva = []
    wstart = p[1:]
    for f in [1.3, PHI, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 10, 12]:
        fun = lambda q: cost(np.concatenate([[np.log(f)], q]), sc)
        starts = [wstart, np.zeros(NW)]
        bestf = min((minimize(fun, s0, method="Nelder-Mead",
                              options=dict(maxfev=160, xatol=1e-3, fatol=1e-6))
                     for s0 in starts), key=lambda o: o.fun)
        curva.append(dict(f=f, c=float(bestf.fun),
                          w=np.exp(np.clip(bestf.x, LB, UB)).tolist()))
        print(f"   L/D={f:5.2f}  migliore={bestf.fun:.5f}", flush=True)

    # forme di confronto
    conf = {}
    for lab, (xx, rr) in {"ellissoide aureo (L/D 1,618)": axi.ellipsoid(PHI),
                          "goccia aurea (L/D 1,618)": axi.body([PHI, 1, 1, 1, 1]),
                          "sfera": axi.ellipsoid(1.0),
                          "goccia standard L/D 4,5": axi.body([4.5, 1, 1, 1, 1])}.items():
        c = cost_shape(xx, rr, sc)
        Dre = sc["Re"] if sc["obiettivo"] == "volume" else sc["Re"] / (2 * rr.max())
        d = axi.drag(xx, rr, Dre, sc["transition"], detail=True)
        conf[lab] = dict(c=float(c), staccato=d["x_sep_turb"] is not None,
                         x_sep=d["x_sep_turb"])
        print(f"   confronto {lab:30s} {c:.5f}  {'(flusso staccato: stima)' if conf[lab]['staccato'] else ''}",
              flush=True)

    return dict(scenario=name, **{k: v for k, v in sc.items()},
                ottimo=dict(c=best, p=p.tolist(), f=det["f"], x_trans=det["x_trans"],
                            x_sep=det["x_sep_turb"], x_spess_max=float(x[i_max] / x[-1]),
                            x=x.tolist(), r=r.tolist()),
                curva=curva, confronto=conf)


if __name__ == "__main__":
    which = sys.argv[1:] or list(SCENARI)
    out = {}
    for n in which:
        out[n] = run(n, SCENARI[n])
        json.dump(out, open("out/axi_risultati.json", "w"), indent=1)
    print("FATTO", flush=True)
