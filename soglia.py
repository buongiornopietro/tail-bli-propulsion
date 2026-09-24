"""Mappa della soglia di convenienza del propulsore in coda (BLI / Goldschmied).

Per ogni (Re_V, regime):
  1. corpo convenzionale ottimo  -> C_conv = C_DV / eta
  2. corpo + fessura + portata ottimi per K = 0, 0.05, 0.10, 0.20
     (K = perdita di pressione totale nel condotto in unità di ½ρU²)
  3. K* = perdita alla quale BLI e convenzionale si equivalgono
Stesso rendimento delle palette (ETA) per elica e pompa: la differenza è
solo fisica del BLI + perdite nel condotto.
"""
import json, sys, time, os
import numpy as np
from multiprocessing import Pool
from scipy.optimize import differential_evolution, minimize

ETA = 0.85
PENDENZA_MAX = None   # limite su |dr/dx| (None = nessun vincolo)
KS = [0.0, 0.05, 0.10, 0.20]
RES = [3e5, 1e6, 3e6, 1e7, 3e7, 1e8]
REGIMI = ["naturale", "forzata"]


def caso(args):
    Re, reg = args[:2]
    seed = args[2] if len(args) > 2 else 0          # 0 = impostazioni originali
    iters = args[3] if len(args) > 3 else 20
    import axi, bli, axi_opt, bli_opt
    import axi_opt as _ao
    _ao.PENDENZA_MAX = PENDENZA_MAX
    axi.CRITERIO = "eN"          # criterio di transizione e^N (Drela-Giles)
    axi.REX_MAX = 1e12           # nessun tetto artificiale: decide e^N
    tag = f"Re{Re:.0e}_{reg}" + (f"_seme{seed}" if seed else "") + (f"_{args[4]}" if len(args) > 4 else "")
    fn = f"out/soglia_{tag}.json"
    if os.path.exists(fn):
        return json.load(open(fn))
    t0 = time.time()
    sc = dict(Re=Re, transition=reg, obiettivo="volume")
    # 1) convenzionale
    rc = differential_evolution(axi_opt.cost, axi_opt.BOUNDS, args=(sc,), popsize=8,
                                maxiter=25 + (iters - 20), tol=1e-4, seed=7 + seed, polish=False, init="sobol")
    lc = minimize(axi_opt.cost, rc.x, args=(sc,), method="Nelder-Mead",
                  options=dict(maxfev=200, xatol=1e-3, fatol=1e-7))
    pc = lc.x if lc.fun < rc.fun else rc.x
    c_conv = min(lc.fun, rc.fun) / ETA
    # 2) BLI per ogni K (partenza calda dall'ottimo del K precedente)
    sc2 = dict(Re=Re, transition=reg)
    out = dict(Re=Re, regime=reg, C_conv=float(c_conv), p_conv=np.asarray(pc).tolist(), bli=[])
    prev = None
    rng = np.random.default_rng(3 + seed)
    for K in KS:
        var = dict(eta_pompa=ETA, K_condotto=K)
        npop = 6 * len(bli_opt.BOUNDS)
        init = bli_opt.LO + rng.random((npop, len(bli_opt.BOUNDS))) * (bli_opt.HI - bli_opt.LO)
        # semina: forma convenzionale ottima + fessura in coda; e l'ottimo precedente
        n1 = 1 + bli_opt.NW
        init[0] = np.concatenate([np.clip(pc, bli_opt.LO[:n1], bli_opt.HI[:n1]), [0.9, np.log(0.06)]])
        if prev is not None:
            init[1] = prev
        r = differential_evolution(bli_opt.cost, bli_opt.BOUNDS, args=(sc2, var), maxiter=iters,
                                   tol=1e-4, seed=11 + seed, polish=False, init=init)
        l = minimize(bli_opt.cost, r.x, args=(sc2, var), method="Nelder-Mead",
                     options=dict(maxfev=200, xatol=1e-3, fatol=1e-7))
        p = np.clip(l.x if l.fun < r.fun else r.x, bli_opt.LO, bli_opt.HI)
        prev = p
        x, rr, xs, Q = bli_opt.unpack(p)
        d = bli.power(x, rr, Re, reg, xs, Q, detail=True, **var)
        out["bli"].append(dict(K=K, CP=float(d["CP"]), guadagno=float(100 * (d["CP"] / c_conv - 1)),
                               f=float(x[-1] / (2 * rr.max())), xs=xs, Q=Q, Q_BL=d["Q_BL"],
                               x_trans=d["x_trans"], pen=d["penalita"], p=p.tolist()))
        print(f"{tag} K={K:.2f}: BLI {d['CP']:.5f} vs conv {c_conv:.5f} ({100 * (d['CP'] / c_conv - 1):+.1f}%) "
              f"L/D {x[-1] / (2 * rr.max()):.1f} fessura {xs:.2f} Q/Q_BL {Q / max(d['Q_BL'], 1e-9):.1f}", flush=True)
    # 3) soglia K*
    Ks = np.array(KS); dlt = np.array([b["CP"] for b in out["bli"]]) - c_conv
    Kstar = None
    if dlt[0] >= 0:
        Kstar = 0.0
    else:
        for i in range(len(Ks) - 1):
            if dlt[i] < 0 <= dlt[i + 1]:
                Kstar = float(Ks[i] + (Ks[i + 1] - Ks[i]) * (-dlt[i]) / (dlt[i + 1] - dlt[i]))
                break
        if Kstar is None:
            Kstar = float("inf")
    out["K_star"] = Kstar
    out["secondi"] = round(time.time() - t0)
    json.dump(out, open(fn, "w"), indent=1)
    print(f"*** {tag}: soglia K* = {Kstar}  ({out['secondi']}s)", flush=True)
    return out


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "ripeti":
        # ripetibilità: stesso caso con semi diversi e più iterazioni
        Re, reg = float(sys.argv[2]), sys.argv[3]
        casi = [(Re, reg, s, 35) for s in (1, 2, 3)]
        with Pool(len(casi)) as pool:
            res = pool.map(caso, casi)
        print("RIPETIZIONI:", [round(100 * r["K_star"], 2) for r in res], flush=True)
        sys.exit()
    casi = [(Re, reg) for reg in REGIMI for Re in RES]
    with Pool(8) as pool:
        res = list(pool.imap_unordered(caso, casi))
    json.dump(res, open("out/soglia_tutti.json", "w"), indent=1)
    print("FATTO", flush=True)
