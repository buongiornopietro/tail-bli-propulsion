"""Griglia strutturata assialsimmetrica (cuneo) con blockMesh per un corpo di
rivoluzione di lunghezza L = 1 con naso in (0,0) e coda in (1,0).

Topologia a C (3 blocchi, piano meridiano x-r):
  blocco A: avancorpo, dal naso al punto di massimo spessore M; il lato
            sinistro è l'asse davanti al naso (collassato)
  blocco B: retrocorpo, da M alla coda
  blocco C: scia, fondo = asse (collassato), dalla coda all'uscita
In direzione normale (j) le linee vanno dal corpo all'arco ellittico esterno.
Il corpo è spezzato in M perché naso e coda stanno sull'asse: i profili del
lato anteriore e posteriore del cuneo avrebbero altrimenti gli stessi estremi.
Il cuneo è ±ALFA gradi attorno all'asse x.
"""
import numpy as np

ALFA = 2.5


def espansione(h0, lung, n):
    """Rapporto r della serie geometrica con primo passo h0, n celle, lunghezza lung."""
    if h0 * n >= lung:
        return 1.0
    lo, hi = 1.0 + 1e-9, 2.0
    while h0 * (hi**n - 1) / (hi - 1) < lung:
        hi *= 1.5
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if h0 * (m**n - 1) / (m - 1) < lung:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


def griglia_parete(h0, spessore_bl, R, n_bl, n_est):
    """Multi-grading normale: n_bl celle fino a spessore_bl con primo passo h0,
    poi n_est celle fino a R che proseguono la crescita."""
    r1 = espansione(h0, spessore_bl, n_bl)
    h_last = h0 * r1 ** (n_bl - 1)
    r2 = espansione(h_last * r1, R - spessore_bl, n_est)
    f1 = spessore_bl / R
    n = n_bl + n_est
    return (f"(({f1:.6f} {n_bl / n:.6f} {r1 ** (n_bl - 1):.6g}) "
            f"({1 - f1:.6f} {n_est / n:.6f} {r2 ** (n_est - 1):.6g}))"), n


def blockMeshDict(x, r, Re_L, yplus=1.0, xin=-4.0, xout=8.0, R=4.0,
                  n_a=120, n_b=160, n_scia=160, n_bl=60, n_est=50):
    Cf = 0.074 * Re_L ** -0.2
    h0 = yplus / (Re_L * np.sqrt(Cf / 2))
    delta = max(0.05, 3 * 0.37 * Re_L ** -0.2)       # zona fitta: ~3 spessori di strato limite
    gj, nj = griglia_parete(h0, delta, R, n_bl, n_est)
    a = np.radians(ALFA)
    im = int(np.argmax(r))
    xm, rm = float(x[im]), float(r[im])
    c = (xm - 1) / (1 - xin)
    thm = np.arccos(np.clip(c, -1, 1))
    tm = (1 + (1 - xin) * np.cos(thm), R * np.sin(thm))

    def v3(xx, rr, lato):          # lato -1 = fronte, +1 = retro
        return f"({xx:.9g} {rr * np.cos(a):.9g} {lato * rr * np.sin(a):.9g})"

    # vertici 2D: 0 naso, 1 coda, 2 (1,R), 3 (xin,0), 4 (xout,0), 5 (xout,R), 6 M, 7 Tm
    P = [(0.0, 0.0), (1.0, 0.0), (1.0, R), (xin, 0.0), (xout, 0.0), (xout, R), (xm, rm), tm]
    sull_asse = [True, True, False, True, True, False, False, False]
    V, idx = [], {}
    for k, (p, ax) in enumerate(zip(P, sull_asse)):
        idx[(k, -1)] = len(V); V.append(v3(*p, -1))
        if ax:
            idx[(k, 1)] = idx[(k, -1)]                 # asse: vertice collassato
        else:
            idx[(k, 1)] = len(V); V.append(v3(*p, 1))
    f = lambda k: idx[(k, -1)]; b = lambda k: idx[(k, 1)]

    def ellisse(t0, t1, n=40):
        th = np.linspace(t0, t1, n)[1:-1]
        return [(1 + (1 - xin) * np.cos(t), R * np.sin(t)) for t in th]

    def poly(k1, k2, pts, lato):
        s = " ".join(v3(px, pr, lato) for px, pr in pts)
        g = f if lato < 0 else b
        return f"    polyLine {g(k1)} {g(k2)} ({s})"

    avan = list(zip(x[1:im], r[1:im]))
    retro = list(zip(x[im + 1:-1], r[im + 1:-1]))
    # lungo il corpo: fitto al naso (blocco A) e alla coda (blocco B)
    gA = "((0.3 0.4 6) (0.7 0.6 3))"         # dal naso: celle che crescono verso M
    gB = "((0.6 0.5 0.4) (0.4 0.5 0.12))"    # da M alla coda: celle che si stringono
    edges_j = f"{gj} {gj} {gj} {gj}"
    txt = "convertToMeters 1;\n\nvertices\n(\n" + "\n".join("    " + v for v in V) + "\n);\n\n"
    txt += f"""blocks
(
    hex ({f(0)} {f(6)} {f(7)} {f(3)} {b(0)} {b(6)} {b(7)} {b(3)}) ({n_a} {nj} 1)
    edgeGrading ({gA} {gA} {gA} {gA}  {edges_j}  1 1 1 1)
    hex ({f(6)} {f(1)} {f(2)} {f(7)} {b(6)} {b(1)} {b(2)} {b(7)}) ({n_b} {nj} 1)
    edgeGrading ({gB} {gB} {gB} {gB}  {edges_j}  1 1 1 1)
    hex ({f(1)} {f(4)} {f(5)} {f(2)} {b(1)} {b(4)} {b(5)} {b(2)}) ({n_scia} {nj} 1)
    edgeGrading (40 40 40 40  {edges_j}  1 1 1 1)
);

edges
(
"""
    for lato in (-1, 1):
        txt += poly(0, 6, avan, lato) + "\n"
        txt += poly(6, 1, retro, lato) + "\n"
        txt += poly(3, 7, ellisse(np.pi, thm), lato) + "\n"
        txt += poly(7, 2, ellisse(thm, np.pi / 2), lato) + "\n"
    txt += ");\n\n"
    txt += "defaultPatch { name defaultFaces; type empty; }\n\n"
    txt += f"""boundary
(
    corpo    {{ type wall;  faces (({f(0)} {f(6)} {b(6)} {b(0)}) ({f(6)} {f(1)} {b(1)} {b(6)})); }}
    esterno  {{ type patch; faces (({f(3)} {f(7)} {b(7)} {b(3)}) ({f(7)} {f(2)} {b(2)} {b(7)})
                                  ({f(2)} {f(5)} {b(5)} {b(2)})); }}
    uscita   {{ type patch; faces (({f(4)} {f(5)} {b(5)} {b(4)})); }}
    front    {{ type wedge; faces (({f(0)} {f(6)} {f(7)} {f(3)}) ({f(6)} {f(1)} {f(2)} {f(7)})
                                  ({f(1)} {f(4)} {f(5)} {f(2)})); }}
    back     {{ type wedge; faces (({b(0)} {b(3)} {b(7)} {b(6)}) ({b(6)} {b(7)} {b(2)} {b(1)})
                                  ({b(1)} {b(2)} {b(5)} {b(4)})); }}
);
"""
    return txt, h0, delta
