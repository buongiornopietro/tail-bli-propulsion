"""Validazione: corpo di Goldschmied (1966), L/D = 3, turbolento, Re_L = 1e7.
Goldschmied: C_D,V totale equivalente 0.0162 (scia 0.0020 + aspirazione 0.0142);
miglior corpo convenzionale (modello Akron) 0.0235."""
import json
import warnings
import numpy as np
from scipy.optimize import differential_evolution
import axi, bli, bli_opt, axi_opt

warnings.filterwarnings("ignore")
ReL = 1e7


def main():
    x0, r0 = axi.body([3.0, 1, 1, 1, 1])
    ReV = ReL / x0[-1]                   # V = 1  ->  Re_V = Re_L V^(1/3) / L
    print(f"L/D=3: L={x0[-1]:.3f} (V=1)  ->  Re_V = {ReV:.3e}", flush=True)

    # BLI ideale (pompa perfetta, condotto senza perdite), L/D fissato a 3
    B = list(bli_opt.BOUNDS)
    B[0] = (np.log(3.0), np.log(3.0) + 1e-9)
    sc = dict(Re=ReV, transition="forzata")
    var = dict(eta_pompa=1.0, K_condotto=0.0)
    r = differential_evolution(bli_opt.cost, B, args=(sc, var), popsize=10, maxiter=30,
                               tol=1e-5, seed=2, polish=False, init="sobol", workers=4,
                               updating="deferred")
    p = np.clip(r.x, bli_opt.LO, bli_opt.HI)
    x, rr, xs, Q = bli_opt.unpack(p)
    d = bli.power(x, rr, ReV, "forzata", xs, Q, detail=True)

    # convenzionale ottimo con L/D = 3, stessa Re_L
    LB = [(np.log(3.0), np.log(3.0) + 1e-9)] + [(axi_opt.LB, axi_opt.UB)] * 4
    rc = differential_evolution(axi_opt.cost, LB,
                                args=(dict(Re=ReV, transition="forzata", obiettivo="volume"),),
                                popsize=10, maxiter=30, seed=2, polish=False, workers=4,
                                updating="deferred")
    out = dict(ReV=ReV, CP_bli=d["CP"], CP_ingestione=d["CP_ingestione"], CP_coda=d["CP_coda"],
               xs=xs, Q=Q, Q_BL=d["Q_BL"], sep_coda=d["sep_coda"], penalita=d["penalita"],
               CDV_conv_LD3=float(rc.fun),
               goldschmied=dict(totale=0.0162, scia=0.0020, aspirazione=0.0142, akron=0.0235),
               x=x.tolist(), r=rr.tolist())
    json.dump(out, open("out/valida_goldschmied.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k not in ("x", "r")}, indent=1), flush=True)


if __name__ == "__main__":
    main()
