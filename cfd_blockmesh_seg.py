"""Griglia strutturata assialsimmetrica (cuneo) con parete del corpo divisa in
segmenti, per rappresentare presa (aspirazione) e scarico (soffiaggio) di un
propulsore in coda "alla Goldschmied".

Topologia a C: una fila di blocchi lungo il corpo (dal naso alla coda), uno per
segmento, più un blocco di scia. Le linee j vanno dalla parete all'arco
ellittico esterno. Naso e coda stanno sull'asse (vertici collassati); ogni
segmento ha almeno un estremo fuori asse, quindi i profili curvi del lato
anteriore e posteriore del cuneo sono distinti.

segmenti: lista di (x_fine/L, n_celle, grading, patch); il primo parte dal naso,
l'ultimo deve finire in 1.0 (coda).
"""
import numpy as np
from cfd_blockmesh import griglia_parete, ALFA


def blockMeshDict(x, r, Re_L, segmenti, yplus=1.0, xin=-4.0, xout=8.0, R=4.0,
                  n_scia=160, n_bl=60, n_est=50):
    Cf = 0.074 * Re_L ** -0.2
    h0 = yplus / (Re_L * np.sqrt(Cf / 2))
    delta = max(0.05, 3 * 0.37 * Re_L ** -0.2)
    gj, nj = griglia_parete(h0, delta, R, n_bl, n_est)
    a = np.radians(ALFA)
    L = x[-1]; x = x / L; r = r / L

    def v3(xx, rr, lato):
        return f"({xx:.9g} {rr * np.cos(a):.9g} {lato * rr * np.sin(a):.9g})"

    def theta(xx):                      # angolo sull'ellisse esterna per l'ascissa xx
        return np.arccos(np.clip((xx - 1) / (1 - xin), -1, 1))

    def arco_pt(t):
        return (1 + (1 - xin) * np.cos(t), R * np.sin(t))

    xs_k = [0.0] + [s[0] for s in segmenti]           # divisioni lungo il corpo
    assert abs(xs_k[-1] - 1.0) < 1e-12
    # vertici 2D: corpo (naso..coda), arco (monte..(1,R)), scia
    P, asse = [], []
    corpo_id, arco_id = [], []
    for k, xx in enumerate(xs_k):
        rr = 0.0 if k in (0, len(xs_k) - 1) else float(np.interp(xx, x, r))
        corpo_id.append(len(P)); P.append((xx, rr)); asse.append(rr == 0.0)
    for k, xx in enumerate(xs_k):
        t = np.pi if k == 0 else (np.pi / 2 if k == len(xs_k) - 1 else theta(xx))
        pt = (xin, 0.0) if k == 0 else arco_pt(t)
        arco_id.append(len(P)); P.append(pt); asse.append(k == 0)
    iw0 = len(P); P.append((xout, 0.0)); asse.append(True)
    iw1 = len(P); P.append((xout, R)); asse.append(False)

    V, idx = [], {}
    for k, (p, ax) in enumerate(zip(P, asse)):
        idx[(k, -1)] = len(V); V.append(v3(*p, -1))
        if ax:
            idx[(k, 1)] = idx[(k, -1)]
        else:
            idx[(k, 1)] = len(V); V.append(v3(*p, 1))
    f = lambda k: idx[(k, -1)]; b = lambda k: idx[(k, 1)]

    ej = f"{gj} {gj} {gj} {gj}"
    blocchi, edges, patch = [], [], {}
    for s, (xe, n, g, nome) in enumerate(segmenti):
        c0, c1 = corpo_id[s], corpo_id[s + 1]
        a0, a1 = arco_id[s], arco_id[s + 1]
        blocchi.append(f"    hex ({f(c0)} {f(c1)} {f(a1)} {f(a0)} {b(c0)} {b(c1)} {b(a1)} {b(a0)}) "
                       f"({n} {nj} 1)\n    edgeGrading ({g} {g} {g} {g}  {ej}  1 1 1 1)")
        x0 = xs_k[s]
        m = (x > x0) & (x < xe)
        pts_c = list(zip(x[m], r[m]))
        t0 = np.pi if s == 0 else theta(x0)
        t1 = np.pi / 2 if s == len(segmenti) - 1 else theta(xe)
        th = np.linspace(t0, t1, 30)[1:-1]
        pts_a = [arco_pt(t) for t in th]
        for lato in (-1, 1):
            gg = f if lato < 0 else b
            if pts_c:
                edges.append(f"    polyLine {gg(c0)} {gg(c1)} (" + " ".join(v3(px, pr, lato) for px, pr in pts_c) + ")")
            edges.append(f"    polyLine {gg(a0)} {gg(a1)} (" + " ".join(v3(px, pr, lato) for px, pr in pts_a) + ")")
        patch.setdefault(nome, []).append(f"({f(c0)} {f(c1)} {b(c1)} {b(c0)})")
        patch.setdefault("esterno", []).append(f"({f(a0)} {f(a1)} {b(a1)} {b(a0)})")
        patch.setdefault("front", []).append(f"({f(c0)} {f(c1)} {f(a1)} {f(a0)})")
        patch.setdefault("back", []).append(f"({b(c0)} {b(a0)} {b(a1)} {b(c1)})")
    # scia
    ct, at = corpo_id[-1], arco_id[-1]
    blocchi.append(f"    hex ({f(ct)} {f(iw0)} {f(iw1)} {f(at)} {b(ct)} {b(iw0)} {b(iw1)} {b(at)}) "
                   f"({n_scia} {nj} 1)\n    edgeGrading (40 40 40 40  {ej}  1 1 1 1)")
    patch["esterno"].append(f"({f(at)} {f(iw1)} {b(iw1)} {b(at)})")
    patch["front"].append(f"({f(ct)} {f(iw0)} {f(iw1)} {f(at)})")
    patch["back"].append(f"({b(ct)} {b(at)} {b(iw1)} {b(iw0)})")
    uscita = [f"({f(iw0)} {f(iw1)} {b(iw1)} {b(iw0)})"]

    tipi = {"front": "wedge", "back": "wedge", "esterno": "patch"}
    txt = "convertToMeters 1;\n\nvertices\n(\n" + "\n".join("    " + v for v in V) + "\n);\n\n"
    txt += "blocks\n(\n" + "\n".join(blocchi) + "\n);\n\n"
    txt += "edges\n(\n" + "\n".join(edges) + "\n);\n\n"
    txt += "defaultPatch { name defaultFaces; type empty; }\n\n"
    txt += "boundary\n(\n"
    for nome, facce in patch.items():
        tp = tipi.get(nome, "wall" if nome == "corpo" else "patch")
        txt += f"    {nome} {{ type {tp}; faces ({' '.join(facce)}); }}\n"
    txt += f"    uscita {{ type patch; faces ({' '.join(uscita)}); }}\n);\n"
    return txt, h0
