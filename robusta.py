"""Forma "garantita": corpo a volume fisso che massimizza il guadagno di velocità MINIMO
rispetto alla goccia turbolenta classica, su un insieme di scenari di transizione.

Perché: la transizione naturale (flusso laminare lungo) può ridurre molto la resistenza,
ma la sua previsione è incerta (criterio, turbolenza esterna, rugosità). Invece di
ottimizzare per uno scenario solo, si cerca la forma che resta migliore in tutti.

Scenari (stesso corpo, stesso Re_V):
  quieto    e^N, Tu = 0,05 %  (N ~ 9,8), tetto Re_x = 1e7 invece di 5e6
  nominale  e^N, Tu = 0,1 %   (N ~ 8,1)
  mosso     e^N, Tu = 0,3 %   (N ~ 5,5)
  agitato   e^N, Tu = 1 %     (N ~ 2,6)
  michel    criterio di Michel (senza storia del gradiente)
  rugoso    e^N, Tu = 0,1 %, transizione comunque a Re_x = 2e6 (rugosità/disturbi)
Riferimento: la forma ottima con strato limite turbolento dal naso (goccia classica),
valutata negli stessi scenari. A potenza fissata la velocità scala come C_DV^(-1/3).

Forme calcolate per ogni Re_V:
  turbolenta  ottimo con transizione forzata (riferimento)
  nominale    ottimo per il solo scenario nominale (quello che si fa di solito)
  robusta     massimizza min_s (C_rif,s / C_s)  -> guadagno garantito

uso: python robusta.py [Re_V ...]      (default 1e6 3e6 1e7 3e7); max 4 processi
"""
import json, sys, time
import numpy as np
from multiprocessing import Pool
from scipy.optimize import differential_evolution, minimize

SCENARI = {
    "quieto":   dict(criterio="eN", Tu=0.0005, Rex_max=1e7),
    "nominale": dict(criterio="eN", Tu=0.001, Rex_max=5e6),
    "mosso":    dict(criterio="eN", Tu=0.003, Rex_max=5e6),
    "agitato":  dict(criterio="eN", Tu=0.01, Rex_max=5e6),
    "michel":   dict(criterio="michel", Tu=0.001, Rex_max=5e6),
    "rugoso":   dict(criterio="eN", Tu=0.001, Rex_max=2e6),
}
DE = dict(popsize=10, maxiter=40, tol=1e-5, seed=5, polish=False, init="sobol")


def valuta(p, ReV, scenari=SCENARI, detail=False):
    """C_DV del corpo di parametri p (log) in ciascuno scenario; il flusso potenziale
    si calcola una volta sola e si riusa."""
    import axi, axi_opt
    p = np.clip(p, [b[0] for b in axi_opt.BOUNDS], [b[1] for b in axi_opt.BOUNDS])
    x, r = axi.body([np.exp(p[0])] + list(np.exp(p[1:])))
    sv = axi.surface_velocity(x, r)
    orig = axi.surface_velocity
    axi.surface_velocity = lambda *_: sv
    out = {}
    try:
        for nome, s in scenari.items():
            axi.CRITERIO, axi.TU = s["criterio"], s["Tu"]
            if s.get("forzata"):
                d = axi.drag(x, r, ReV, "forzata", detail=True)
            else:
                d = axi.drag(x, r, ReV, "naturale", detail=True, Rex_max=s["Rex_max"])
            out[nome] = d if detail else d["CDV"]
    finally:
        axi.surface_velocity = orig
        axi.CRITERIO, axi.TU = "eN", 0.001
    return out


def _ottimizza(f, args):
    import axi_opt
    res = differential_evolution(f, axi_opt.BOUNDS, args=args, **DE)
    loc = minimize(f, res.x, args=args, method="Nelder-Mead",
                   options=dict(maxfev=400, xatol=1e-4, fatol=1e-7))
    return (loc.x, loc.fun) if loc.fun < res.fun else (res.x, res.fun)


def c_turb(p, ReV):
    return valuta(p, ReV, {"t": dict(criterio="eN", Tu=0.001, Rex_max=5e6, forzata=True)})["t"]


def c_nom(p, ReV):
    return valuta(p, ReV, {"n": SCENARI["nominale"]})["n"]


def c_rob(p, ReV, C_rif):
    c = valuta(p, ReV)
    return max(c[s] / C_rif[s] for s in SCENARI)        # peggior rapporto = da minimizzare


def caso(ReV):
    import warnings; warnings.filterwarnings("ignore")
    t0 = time.time()
    p_t, _ = _ottimizza(c_turb, (ReV,))
    C_rif = valuta(p_t, ReV)
    p_n, _ = _ottimizza(c_nom, (ReV,))
    p_r, _ = _ottimizza(c_rob, (ReV, C_rif))
    out = dict(ReV=ReV, secondi=time.time() - t0, C_rif=C_rif, forme={})
    for nome, p in (("turbolenta", p_t), ("nominale", p_n), ("robusta", p_r)):
        det = valuta(p, ReV, detail=True)
        import axi, axi_opt
        pc = np.clip(p, [b[0] for b in axi_opt.BOUNDS], [b[1] for b in axi_opt.BOUNDS])
        x, r = axi.body([np.exp(pc[0])] + list(np.exp(pc[1:])))
        out["forme"][nome] = dict(
            p=list(map(float, p)), LsuD=det["nominale"]["f"],
            x_spess_max=float(x[np.argmax(r)] / x[-1]),
            C={s: det[s]["CDV"] for s in SCENARI},
            x_trans={s: det[s]["x_trans"] for s in SCENARI},
            sep_lam={s: det[s]["sep_lam"] for s in SCENARI},
            x_sep_turb={s: det[s]["x_sep_turb"] for s in SCENARI},
            guadagno_vel={s: 100 * ((C_rif[s] / det[s]["CDV"]) ** (1 / 3) - 1) for s in SCENARI},
            x=list(map(float, x / x[-1])), r=list(map(float, r / x[-1])))
    return out


if __name__ == "__main__":
    Res = [float(a) for a in sys.argv[1:]] or [1e6, 3e6, 1e7, 3e7]
    with Pool(min(4, len(Res))) as pool:
        ris = list(pool.imap_unordered(caso, Res))
    ris.sort(key=lambda d: d["ReV"])
    for d in ris:
        print(f"\nRe_V = {d['ReV']:.0e}  ({d['secondi']/60:.0f} min)")
        for nome, f in d["forme"].items():
            g = f["guadagno_vel"]
            print(f"  {nome:10s} L/D {f['LsuD']:5.2f}  spess.max a {f['x_spess_max']:.2f}  "
                  "vel. vs goccia: " + "  ".join(f"{s} {g[s]:+5.1f}%" for s in SCENARI)
                  + f"  | minimo {min(g.values()):+.1f}%", flush=True)
    json.dump(ris, open("out/robusta.json", "w"), indent=1)
    print("\nout/robusta.json")
