"""Ottimizzazione automatica della forma a resistenza minima, a PARI AREA.

Forma: profilo simmetrico CST (Kulfan)
    y(ξ) = ξ^0.5 (1-ξ) Σ_k w_k B_k(ξ),   ξ = x/L ∈ [0,1]
(muso arrotondato, coda chiusa). Parametri ottimizzati:
    p[0]   = ln(L/H)        snellezza
    p[1:]  = ln(w_k/w_0)    distribuzione dello spessore (w_0 = 1)
La forma viene poi scalata perché l'area sia sempre A0.

Algoritmo: strategia evolutiva (μ/μ_w, λ) con passo adattivo, λ valutazioni
in parallelo per generazione.

uso: python opt.py [generazioni]
"""
import sys, json, time
import numpy as np
from math import comb
from multiprocessing import Pool
from lbm import simulate, make_shape, spiral_mask

A0 = 600            # area del corpo (celle)  -> "volume" fisso
Re = 160            # Reynolds basato su sqrt(A0): tau = 0.545 (collaudato stabile)
U = 0.1             # Mach basso: con U = 0.12 divergeva
NU = U * np.sqrt(A0) / Re
NX, NY = 500, 150
X0 = 120          # ~5 H dall ingresso: vicino all ingresso diverge
STEPS = int(2.6 * NX / U)   # ~2.6 tempi di attraversamento
AVG = 0.55                  # media delle forze sull'ultimo 45%
SPONGE = 0.0          # la zona spugna destabilizzava il TRT: disattivata
NW = 4                      # numero di pesi di Bernstein
LAM, MU = 8, 4
LOG_F = (np.log(1.3), np.log(9.0))

xi = np.linspace(0, 1, 801)


def profile(p):
    """Semispessore adimensionale y(ξ) (spessore max = 1/(L/H) con L = 1)."""
    fr = np.exp(np.clip(p[0], *LOG_F))
    wts = np.concatenate([[1.0], np.exp(np.clip(p[1:], -3, 3))])
    n = len(wts) - 1
    B = sum(wts[k] * comb(n, k) * xi**k * (1 - xi)**(n - k) for k in range(n + 1))
    y = np.sqrt(xi) * (1 - xi) * B
    return y / (2 * y.max()) / fr, fr


def mask_from_params(p):
    y, fr = profile(p)
    area_unit = 2 * np.trapezoid(y, xi)
    L = np.sqrt(A0 / area_unit)
    Y, X = np.mgrid[0:NY, 0:NX].astype(float) + 0.5
    xs = (X - X0) / L
    ok = (xs >= 0) & (xs <= 1)
    yt = np.interp(np.clip(xs, 0, 1), xi, y) * L
    solid = ok & (np.abs(Y - NY / 2 - 0.3) <= yt)
    return solid, L, L / fr


def evaluate(args):
    tag, solid = args
    t0 = time.time()
    Fx, Fy, _ = simulate(solid, NU, U, STEPS, avg_from=AVG, sponge=SPONGE, outlet="pressione")
    q = 0.5 * U**2
    ref = q * np.sqrt(A0)
    cd = float(Fx.mean() / ref)
    if not np.isfinite(cd):
        cd = 99.0
    return dict(tag=tag, Cd_area=cd, Cl=float(Fy.mean() / ref),
                osc=float(Fx.std() / ref), area=int(solid.sum()),
                sec=round(time.time() - t0))


def reference_masks():
    """Forme di confronto, riscalate alla stessa area A0."""
    Y, X = np.mgrid[0:NY, 0:NX].astype(float) + 0.5
    out = {}
    for name, a_unit in [("cerchio", np.pi / 4), ("ellisse_aurea", np.pi / 4 * 1.618034),
                         ("goccia_aurea", None), ("spirale_aurea", None)]:
        if name == "spirale_aurea":
            m1, _ = spiral_mask(X, Y, X0, NY / 2 + 0.3, 40.0)
            H = 40.0 * np.sqrt(A0 / m1.sum())
            m, _ = spiral_mask(X, Y, X0, NY / 2 + 0.3, H)
        elif name == "goccia_aurea":
            m1, _ = make_shape(name, X, Y, X0, NY / 2 + 0.3, 40.0)
            H = 40.0 * np.sqrt(A0 / m1.sum())
            m, _ = make_shape(name, X, Y, X0, NY / 2 + 0.3, H)
        else:
            H = np.sqrt(A0 / a_unit)
            m, _ = make_shape(name, X, Y, X0, NY / 2 + 0.3, H)
        out[name] = m
    return out


def main(gens):
    rng = np.random.default_rng(1)
    dim = 1 + (NW - 1)
    # partenza: goccia NACA-like con L/H = 3
    mean = np.array([np.log(3.0)] + [0.0] * (NW - 1))
    sigma = np.array([0.35] + [0.6] * (NW - 1))
    wts = np.log(MU + 0.5) - np.log(np.arange(1, MU + 1))
    wts /= wts.sum()

    history = []
    with Pool(LAM) as pool:
        best = None
        for g in range(gens):
            z = rng.standard_normal((LAM, dim))
            z[LAM // 2:] = -z[:LAM // 2]            # campionamento speculare
            cand = mean + sigma * z
            cand[:, 0] = np.clip(cand[:, 0], *LOG_F)
            jobs = []
            for k, p in enumerate(cand):
                solid, L, H = mask_from_params(p)
                jobs.append((f"g{g}_{k}", solid))
            res = pool.map(evaluate, jobs)
            cds = np.array([r["Cd_area"] for r in res])
            order = np.argsort(cds)
            for k, (p, r) in enumerate(zip(cand, res)):
                _, L, H = mask_from_params(p)
                history.append(dict(gen=g, p=p.tolist(), LsuH=float(L / H),
                                    L=float(L), H=float(H), **r))
            if best is None or cds[order[0]] < best["Cd_area"]:
                best = history[-LAM + order[0]]
            old = mean.copy()
            mean = (wts[:, None] * cand[order[:MU]]).sum(0)
            # passo adattivo: se la media migliora molto allarga, altrimenti restringe
            succ = (cds < np.median(cds)).mean()
            sigma *= 0.82 if g > 1 else 0.95
            print(f"gen {g}: migliore Cd_area={cds[order[0]]:.4f} "
                  f"L/H={history[-LAM + order[0]]['LsuH']:.2f}  "
                  f"media L/H={np.exp(mean[0]):.2f}  "
                  f"(assoluto migliore {best['Cd_area']:.4f}, L/H={best['LsuH']:.2f})",
                  flush=True)
            json.dump(dict(history=history, best=best, mean=mean.tolist(),
                           sigma=sigma.tolist()),
                      open("out/opt_storia.json", "w"), indent=1)

        # verifica finale della media della distribuzione (spesso più robusta)
        solid, L, H = mask_from_params(mean)
        fin = pool.map(evaluate, [("media_finale", solid)])[0]
        fin.update(p=mean.tolist(), LsuH=float(L / H), L=float(L), H=float(H))
        print("media finale:", json.dumps(fin), flush=True)
        json.dump(dict(history=history, best=best, mean=mean.tolist(),
                       finale=fin, sigma=sigma.tolist()),
                  open("out/opt_storia.json", "w"), indent=1)
    print("FATTO", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 8)
