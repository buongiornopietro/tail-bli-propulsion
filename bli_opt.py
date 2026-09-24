"""Ottimizzazione congiunta forma + fessura + portata per corpo con propulsore in coda.

Variabili (7): ln(L/D), ln w_1..w_4 (CST), x_s/L, ln Q.
Varianti: 'ideale' (pompa perfetta, condotto senza perdite) e
'realistica' (pompa 85%, perdita condotto K = 0,10·½ρU² sulla portata;
riferimento convenzionale con elica all'80%).
"""
import json, sys, time
import numpy as np
from scipy.optimize import differential_evolution, minimize
import axi, bli

NW = 4
BOUNDS = ([(np.log(1.3), np.log(12.0))] + [(np.log(0.1), np.log(5.0))] * NW
          + [(0.55, 0.96), (np.log(0.01), np.log(0.25))])
LO = np.array([b[0] for b in BOUNDS]); HI = np.array([b[1] for b in BOUNDS])


def set_nw(n):
    """Cambia il numero di pesi CST (flessibilità della forma)."""
    global NW, BOUNDS, LO, HI
    NW = n
    BOUNDS = ([(np.log(1.3), np.log(12.0))] + [(np.log(0.1), np.log(5.0))] * NW
              + [(0.55, 0.96), (np.log(0.01), np.log(0.25))])
    LO = np.array([b[0] for b in BOUNDS]); HI = np.array([b[1] for b in BOUNDS])

SCEN = {
    "A_turbolento_ReV1e7": dict(Re=1e7, transition="forzata"),
    "C_naturale_ReV1e6": dict(Re=1e6, transition="naturale"),
}
VAR = {"ideale": dict(eta_pompa=1.0, K_condotto=0.0, eta_elica=1.0),
       "realistica": dict(eta_pompa=0.85, K_condotto=0.10, eta_elica=0.80)}


N_LOCALE = 60        # pannelli aggiuntivi concentrati attorno alla presa


def unpack(p):
    p = np.clip(p, LO, HI)
    xs = float(p[1 + NW])
    xi = axi.nodes(axi.N_PANNELLI, x_fitto=xs, extra=N_LOCALE)
    x, r = axi.body([np.exp(p[0])] + list(np.exp(p[1:1 + NW])), xi=xi)
    return x, r, xs, float(np.exp(p[2 + NW]))


def cost(p, sc, var):
    try:
        x, r, xs, Q = unpack(p)
        import axi_opt
        pen = axi_opt.penalita_forma(x, r)
        c = bli.power(x, r, sc["Re"], sc["transition"], xs, Q,
                      eta_pompa=var["eta_pompa"], K_condotto=var["K_condotto"])
        return float(c + pen) if np.isfinite(c) else 10.0
    except Exception:
        return 10.0


def run(sn, vn):
    sc, var = SCEN[sn], VAR[vn]
    t0 = time.time()
    res = differential_evolution(cost, BOUNDS, args=(sc, var), popsize=8, maxiter=30,
                                 tol=1e-4, seed=5, polish=False, init="sobol")
    loc = minimize(cost, res.x, args=(sc, var), method="Nelder-Mead",
                   options=dict(maxfev=250, xatol=1e-3, fatol=1e-6))
    p = np.clip(loc.x if loc.fun < res.fun else res.x, LO, HI)
    x, r, xs, Q = unpack(p)
    d = bli.power(x, r, sc["Re"], sc["transition"], xs, Q, detail=True,
                  eta_pompa=var["eta_pompa"], K_condotto=var["K_condotto"])
    # riferimento: miglior corpo convenzionale (da axi_opt) con elica
    ref = json.load(open("out/axi_risultati.json"))[sn]["ottimo"]["c"] / var["eta_elica"]
    # BLI applicato al corpo ottimo convenzionale (solo x_s e Q ottimizzati)
    o = json.load(open("out/axi_risultati.json"))[sn]["ottimo"]
    xo, ro = np.array(o["x"]), np.array(o["r"])
    f2 = lambda q: bli.power(xo, ro, sc["Re"], sc["transition"], np.clip(q[0], .55, .96),
                             np.exp(np.clip(q[1], LO[-1], HI[-1])),
                             eta_pompa=var["eta_pompa"], K_condotto=var["K_condotto"])
    r2 = differential_evolution(f2, [(0.55, 0.96), (LO[-1], HI[-1])], popsize=10,
                                maxiter=25, seed=1, polish=False)
    i_max = int(np.argmax(r))
    out = dict(scenario=sn, variante=vn, CP=d["CP"], rif_convenzionale=ref,
               guadagno=100 * (d["CP"] / ref - 1),
               CP_bli_su_corpo_conv=float(r2.fun), xs_su_corpo_conv=float(r2.x[0]),
               Q_su_corpo_conv=float(np.exp(r2.x[1])),
               f=float(x[-1] / (2 * r.max())), x_spess_max=float(x[i_max] / x[-1]),
               xs=xs, Q=Q, Q_BL=d["Q_BL"], CP_ingestione=d["CP_ingestione"],
               CP_coda=d["CP_coda"], penalita=d["penalita"], sep_avancorpo=d["sep_avancorpo"],
               sep_coda=d["sep_coda"], x_trans=d["x_trans"], Ue_fessura=d["Ue_fessura"],
               p=p.tolist(), x=x.tolist(), r=r.tolist(), secondi=round(time.time() - t0))
    print(f"[{sn} | {vn}] C_P={d['CP']:.5f} vs convenzionale {ref:.5f} ({out['guadagno']:+.1f}%)  "
          f"L/D={out['f']:.2f} spess.max {out['x_spess_max']:.2f}  fessura x/L={xs:.2f} Q={Q:.3f} "
          f"(Q_BL {d['Q_BL']:.3f})  sep coda={d['sep_coda']}  | BLI su corpo conv: {r2.fun:.5f} "
          f"({100 * (r2.fun / ref - 1):+.1f}%)  [{out['secondi']}s]", flush=True)
    return out


if __name__ == "__main__":
    res = {}
    for sn in SCEN:
        for vn in VAR:
            res[f"{sn}|{vn}"] = run(sn, vn)
            json.dump(res, open("out/bli_risultati.json", "w"), indent=1)
    print("FATTO", flush=True)
