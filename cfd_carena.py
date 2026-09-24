"""Propulsore in coda INTUBATO: carenatura anulare con labbro arrotondato attorno
alla coda del corpo, ventilatore (disco attuatore) dentro il condotto.

Geometria (unità di L, corpo dal naso (0,0) alla coda (1,0)):
  carenatura = profilo NACA a 4 cifre simmetrico, spessore relativo TAU, corda C,
  bordo d'attacco in X_LE, linea media rettilinea da R_LE a R_TE (convergente:
  l'uscita fa da ugello). Spessore applicato in verticale (sezione "tranciata").

Griglia blockMesh a 3 strati attorno al corpo (cuneo di 5 gradi):
  strato I   : parete del corpo -> curva C1 (= superficie interna della carenatura)
  strato II  : C1 -> C2 (= superficie esterna); nel tratto della carenatura è vuoto
  strato III : C2 -> arco esterno
  colonne    : naso | max spessore | spigoli del labbro | bordo d'uscita | coda | uscita
Il labbro (tratto del contorno con pendenza > 45°) è la faccia di valle del blocco
II a monte; il bordo d'uscita (tronco, spessore TTE) è la faccia di monte del blocco
II a valle.

uso: python cfd_carena.py <nome> [T_iniziale]
"""
import json, os, re, sys
import numpy as np
import cfd_caso
from cfd_blockmesh import ALFA

CHIAVE = "conv_1e6f"                 # corpo ottimo turbolento Re_V = 1e6
X_LE, C, TAU, TTE = 0.855, 0.095, 0.05, 0.0006
R_LE, R_TE = 0.0293, 0.0249
X_F, SP_F = 0.905, 0.010             # ventilatore: centro e spessore del disco
XIN, XOUT, R = -4.0, 8.0, 4.0
H_CORPO_YP = 0.7                     # y+ del primo strato sul corpo
H_CAR = 1.0e-5                       # primo strato sulla carenatura (y+ ~ 2)
FATTORE = 1.0                        # moltiplica il numero di celle in ogni direzione
G1N, T2N, T2M = 0.0275, 0.02, 0.01    # spessori delle fasce davanti al naso e al max spessore
X_T = 0.70                           # colonna di transizione: da qui in poi fasce fitte


def yt(s):
    """Semispessore NACA 4 cifre (bordo chiuso) + bordo d'uscita tronco TTE."""
    s = np.asarray(s, float)
    return (5 * TAU * C * (0.2969 * np.sqrt(s) - 0.1260 * s - 0.3516 * s**2 + 0.2843 * s**3 - 0.1015 * s**4)
            + 0.5 * TTE * s)


def r_media(x):
    return R_LE + (R_TE - R_LE) * (np.asarray(x) - X_LE) / C


def carena(x):
    s = np.clip((np.asarray(x) - X_LE) / C, 0, 1)
    return r_media(x) - yt(s), r_media(x) + yt(s)


def s_spigolo():
    """Ascissa di corda dove la pendenza del contorno scende a 45 gradi."""
    s = np.geomspace(1e-6, 0.2, 20000)
    d = np.gradient(yt(s), s * C)
    return float(s[np.argmax(d < 1.0)])


def grading(L, n, ha, hb=None):
    """Multi-grading blockMesh: n celle su lunghezza L, prima cella ha, ultima hb
    (se hb è None la crescita è a senso unico)."""
    if hb is None:
        lo, hi = 1.0 + 1e-12, 3.0
        for _ in range(200):
            q = 0.5 * (lo + hi)
            if ha * (q**n - 1) / (q - 1) < L: lo = q
            else: hi = q
        return f"{q ** (n - 1):.6g}"

    def lung(q):
        a = np.log(hb / ha) / np.log(q)
        n1 = np.clip(round((n + a) / 2), 1, n - 1); n2 = n - n1
        l1 = ha * (q**n1 - 1) / (q - 1); l2 = hb * (q**n2 - 1) / (q - 1)
        return l1, l2, n1, n2

    lo, hi = 1.0 + 1e-9, 3.0
    for _ in range(200):
        q = 0.5 * (lo + hi)
        l1, l2, _, _ = lung(q)
        if l1 + l2 < L: lo = q
        else: hi = q
    l1, l2, n1, n2 = lung(q)
    f1 = l1 / (l1 + l2)
    return (f"(({f1:.6f} {n1 / n:.6f} {q ** (n1 - 1):.6g}) "
            f"({1 - f1:.6f} {n2 / n:.6f} {1 / q ** (n2 - 1):.6g}))")


def blockMeshDict(x, r, Re_L):
    Cf = 0.074 * Re_L ** -0.2
    h0 = H_CORPO_YP / (Re_L * np.sqrt(Cf / 2))
    a = np.radians(ALFA)
    sc = s_spigolo(); xc = X_LE + sc * C; xte = X_LE + C
    im = int(np.argmax(r)); xm, rm = float(x[im]), float(r[im])
    rb = lambda xx: float(np.interp(xx, x, r))
    ri_c, ro_c = (float(v) for v in carena(xc))
    ri_t, ro_t = (float(v) for v in carena(xte))
    g1 = ri_c - rb(xc); t2 = ro_c - ri_c                  # altezza condotto e spessore labbro

    def v3(xx, rr, lato):
        return f"({xx:.10g} {rr * np.cos(a):.10g} {lato * rr * np.sin(a):.10g})"

    def theta(xx):
        return np.arccos(np.clip((xx - 1) / (1 - XIN), -1, 1))

    def arco(t):
        return (1 + (1 - XIN) * np.cos(t), R * np.sin(t))

    # curve di divisione a monte: avancorpo = offset normale; retrocorpo = offset verticale
    m0 = x <= xm
    nx_, nr_ = -np.gradient(r[m0], x[m0]), np.ones(m0.sum())
    nn = np.hypot(nx_, nr_); nx_ /= nn; nr_ /= nn
    nx_[0], nr_[0] = -1.0, 0.0                              # naso: normale assiale
    # davanti al naso le fasce sono spesse (G1N, T2N): niente celle degeneri sull'asse;
    # si assottigliano verso la presa, dove C1 e C2 coincidono con la carenatura
    w0 = (1 - x[m0] / xm) ** 2
    d1_0 = g1 + (G1N - g1) * w0; d2_0 = d1_0 + T2M + (T2N - T2M) * w0
    d2 = lambda xx: g1 + t2 + (T2M - t2) * np.clip((X_T - xx) / (X_T - xm), 0, 1) ** 2
    def offs(liv, zona):
        if zona == 0:
            d = d1_0 if liv == 1 else d2_0
            return list(zip(x[m0] + d * nx_, r[m0] + d * nr_))
        lo, hi = (xm, X_T) if zona == 1 else (X_T, xc)
        m = (x > lo) & (x < hi)
        return list(zip(x[m], r[m] + (g1 if liv == 1 else d2(x[m]))))

    # colonne: 0 naso, 1 max spessore, 2 transizione, 3 spigoli labbro, 4 bordo d'uscita,
    # 5 coda, 6 uscita
    X = [0.0, xm, X_T, xc, xte, 1.0]
    liv = {  # (colonna, livello) -> punto 2D; livelli 0 parete/asse, 1 C1, 2 C2, 3 esterno
        (0, 0): (0.0, 0.0), (0, 1): (-G1N, 0.0), (0, 2): (-G1N - T2N, 0.0), (0, 3): (XIN, 0.0),
        (1, 0): (xm, rm), (1, 1): (xm, rm + g1), (1, 2): (xm, rm + g1 + T2M), (1, 3): arco(theta(xm)),
        (2, 0): (X_T, rb(X_T)), (2, 1): (X_T, rb(X_T) + g1), (2, 2): (X_T, rb(X_T) + g1 + t2),
        (2, 3): arco(theta(X_T)),
        (3, 0): (xc, rb(xc)), (3, 1): (xc, ri_c), (3, 2): (xc, ro_c), (3, 3): arco(theta(xc)),
        (4, 0): (xte, rb(xte)), (4, 1): (xte, ri_t), (4, 2): (xte, ro_t), (4, 3): arco(theta(xte)),
        (5, 0): (1.0, 0.0), (5, 1): (1.0, ri_t), (5, 2): (1.0, ro_t), (5, 3): (1.0, R),
        (6, 0): (XOUT, 0.0), (6, 1): (XOUT, ri_t), (6, 2): (XOUT, ro_t), (6, 3): (XOUT, R),
    }
    V, idx = [], {}
    for key, (px, pr) in liv.items():
        idx[(key, -1)] = len(V); V.append(v3(px, pr, -1))
        if pr == 0.0:
            idx[(key, 1)] = idx[(key, -1)]
        else:
            idx[(key, 1)] = len(V); V.append(v3(px, pr, 1))
    F = lambda c, l: idx[((c, l), -1)]; B = lambda c, l: idx[((c, l), 1)]

    # linee curve (i = lungo x, tra colonne c e c+1 al livello l)
    curve = {}
    xx0 = x[(x > 0) & (x < xm)]; curve[(0, 0)] = list(zip(xx0, r[(x > 0) & (x < xm)]))
    for c in (1, 2, 3, 4):
        m = (x > X[c]) & (x < X[c + 1]); curve[(c, 0)] = list(zip(x[m], r[m]))
    curve[(0, 1)] = offs(1, 0)[1:-1]; curve[(1, 1)] = offs(1, 1); curve[(2, 1)] = offs(1, 2)
    curve[(0, 2)] = offs(2, 0)[1:-1]; curve[(1, 2)] = offs(2, 1); curve[(2, 2)] = offs(2, 2)
    xs_c = xc + (xte - xc) * (1 - np.cos(np.linspace(0, np.pi, 120)[1:-1])) / 2
    ri, ro = carena(xs_c)
    curve[(3, 1)] = list(zip(xs_c, ri)); curve[(3, 2)] = list(zip(xs_c, ro))
    for c in range(5):
        t0 = np.pi if c == 0 else theta(X[c]); t1 = np.pi / 2 if c == 4 else theta(X[c + 1])
        curve[(c, 3)] = [arco(t) for t in np.linspace(t0, t1, 40)[1:-1]]
    # labbro: dallo spigolo interno, attorno al bordo d'attacco, allo spigolo esterno
    ss = sc * np.sin(np.linspace(0, np.pi / 2, 50)[1:-1]) ** 2      # fitto al bordo d'attacco
    xl = X_LE + ss * C; yl = yt(ss); rl = r_media(xl)
    labbro = list(zip(xl[::-1], (rl - yl)[::-1])) + [(X_LE, float(r_media(X_LE)))] + \
             list(zip(xl, rl + yl))

    edges = []
    def poly(k1, k2, pts, lato):
        if pts:
            edges.append(f"    polyLine {k1} {k2} (" + " ".join(v3(px, pr, lato) for px, pr in pts) + ")")
    for (c, l), pts in curve.items():
        poly(F(c, l), F(c + 1, l), pts, -1); poly(B(c, l), B(c + 1, l), pts, 1)
    poly(F(3, 1), F(3, 2), labbro, -1); poly(B(3, 1), B(3, 2), labbro, 1)

    # gradazioni
    n = lambda k: int(round(k * FATTORE))             # infittimento della griglia
    colonne = {0: (n(100), grading(0.33, n(100), 3e-4, 4e-3)),
               1: (n(80), grading(X_T - xm, n(80), 4e-3, 3e-3)),
               2: (n(80), grading(xc - X_T, n(80), 3e-3, H_CAR)),
               3: (n(170), grading(xte - xc, n(170), H_CAR, 6e-5)),
               4: (n(60), grading(1 - xte, n(60), 6e-5, 2e-3)),
               5: (n(140), grading(XOUT - 1, n(140), 2e-3))}
    NJ = {0: n(100), 1: n(36), "1v": n(16), 2: n(90)}
    def jspec(cp, l):
        """Gradazione dei lati j sulla colonna cp: rada al naso e al max spessore,
        fitta su entrambe le pareti dalla presa in poi."""
        if cp <= 1:
            G, T = (G1N, T2N) if cp == 0 else (g1, T2M)
            return {0: grading(G, n(100), h0), 1: "1", 2: grading(R - G - T, n(90), T / n(36))}[l]
        if l == 1 and cp >= 4:
            return grading(TTE, n(16), 1.5e-5, 1.5e-5)
        return {0: grading(g1, n(100), h0, H_CAR), 1: grading(t2, n(36), H_CAR, H_CAR),
                2: grading(R - ro_c, n(90), H_CAR)}[l]
    blocchi, patch = [], {}
    def add(nome, facce): patch.setdefault(nome, []).append(facce)
    for c in range(6):
        for l in range(3):
            if c == 3 and l == 1:
                continue                                   # dentro la carenatura
            rk = "1v" if (l == 1 and c >= 4) else l
            ni, gi = colonne[c]; nj = NJ[rk]
            gl, gr = jspec(c, l), jspec(c + 1, l)
            v = [F(c, l), F(c + 1, l), F(c + 1, l + 1), F(c, l + 1),
                 B(c, l), B(c + 1, l), B(c + 1, l + 1), B(c, l + 1)]
            # lato esterno (arco) a passo uniforme: le celle fitte vicino al labbro
            # non devono propagarsi fino al contorno lontano
            ga = "1" if l == 2 else gi
            blocchi.append(f"    hex ({' '.join(map(str, v))}) ({ni} {nj} 1)\n"
                           f"    edgeGrading ({gi} {ga} {ga} {gi}  {gl} {gr} {gr} {gl}  1 1 1 1)")
            add("front", f"({v[0]} {v[1]} {v[2]} {v[3]})")
            add("back", f"({v[4]} {v[7]} {v[6]} {v[5]})")
            if l == 0 and c < 5:
                add("corpo", f"({v[0]} {v[1]} {v[5]} {v[4]})")
            if l == 2:
                add("esterno", f"({v[3]} {v[2]} {v[6]} {v[7]})")
            if c == 5:
                add("uscita", f"({v[1]} {v[2]} {v[6]} {v[5]})")
            if (c == 3 and l == 0):
                add("carena", f"({v[3]} {v[2]} {v[6]} {v[7]})")     # superficie interna
            if (c == 3 and l == 2):
                add("carena", f"({v[0]} {v[1]} {v[5]} {v[4]})")     # superficie esterna
            if (c == 2 and l == 1):
                add("carena", f"({v[1]} {v[2]} {v[6]} {v[5]})")     # labbro
            if (c == 4 and l == 1):
                add("carena", f"({v[0]} {v[3]} {v[7]} {v[4]})")     # bordo d'uscita
    tipi = {"front": "wedge", "back": "wedge", "corpo": "wall", "carena": "wall"}
    txt = "convertToMeters 1;\n\nvertices\n(\n" + "\n".join("    " + v for v in V) + "\n);\n\n"
    txt += "blocks\n(\n" + "\n".join(blocchi) + "\n);\n\n"
    txt += "edges\n(\n" + "\n".join(edges) + "\n);\n\n"
    txt += "defaultPatch { name defaultFaces; type empty; }\n\n"
    txt += "boundary\n(\n"
    for nome, facce in patch.items():
        txt += f"    {nome} {{ type {tipi.get(nome, 'patch')}; faces ({' '.join(facce)}); }}\n"
    txt += ");\n"
    geo = dict(xc=xc, xte=xte, s_spigolo=sc, g1=g1, spessore_labbro=t2, ri_c=ri_c, ro_c=ro_c,
               ri_te=ri_t, ro_te=ro_t, h0=h0, raggio_LE=1.1019 * TAU**2 * C)
    return txt, geo


def main(nome, T0=7.2e-6):
    prof = json.load(open("out/cfd_profili.json"))[CHIAVE]
    x, r = np.array(prof["x"]), np.array(prof["r"]); L = x[-1]; x = x / L; r = r / L
    Re_L = prof["ReL"]
    cfd_caso.caso(nome, x, r, Re_L, H_CORPO_YP, modello="kOmegaSST")
    bm, geo = blockMeshDict(x, r, Re_L)
    cfd_caso.scrivi(nome, "system/blockMeshDict", "dictionary", bm)
    # condizioni al contorno della carenatura = parete (funzione di parete di Spalding)
    for campo in ("U", "p", "k", "omega", "nut"):
        fn = os.path.join(nome, "0", campo)
        t = open(fn).read()
        riga = re.search(r"\n(\s*)corpo(\s*\{[^\n]*\})", t)
        bc = riga.group(2)
        if campo == "nut":
            bc = " { type nutUSpaldingWallFunction; value uniform 0; }"
        t = t[:riga.end()] + f"\n{riga.group(1)}carena{bc}" + t[riga.end():]
        open(fn, "w", newline="\n").write(t)
    # disco attuatore nel condotto: cilindro di raggio pari alla linea media
    rf = float(r_media(X_F))
    zona = f"""cellZone
        {{
            type            cylinder;
            point1          ({X_F - SP_F / 2} 0 0);
            point2          ({X_F + SP_F / 2} 0 0);
            radius          {rf:.6f};
        }}"""
    cfd_caso.scrivi(nome, "constant/fvModels", "dictionary", f"""
    disco
    {{
        type            semiImplicitSource;
        {zona}
        volumeMode      absolute;
        sources {{ U {{ explicit ({T0:.6e} 0 0); implicit 0; }} }}
    }}
    """)
    cd = open(os.path.join(nome, "system/controlDict")).read()
    cd = re.sub(r"patches(\s+)\(corpo\);(\s+rho)", r"patches\1(corpo carena);\2", cd, count=1)
    extra = f"""
        forze_carena
        {{
            type            forces;
            libs            ("libforces.so");
            writeControl    timeStep;
            writeInterval   10;
            patches         (carena);
            rho             rhoInf;
            rhoInf          1;
            CofR            (0 0 0);
        }}
        Udisco
        {{
            type            volFieldValue;
            libs            ("libfieldFunctionObjects.so");
            writeControl    timeStep;
            writeInterval   10;
            writeFields     false;
            {zona}
            operation       volAverage;
            fields          (U);
        }}
        ptot
        {{
            type            pressure;
            libs            ("libfieldFunctionObjects.so");
            calcTotal       yes;
            calcCoeff       no;
            rho             rhoInf;
            rhoInf          1;
            pRef            0;
            executeControl  writeTime;
            writeControl    writeTime;
        }}
"""
    i = re.search(r"\n\s*yPlus \{", cd).start()
    cd = cd[:i] + "\n" + extra + cd[i:]
    cd = cd.replace("patches (corpo); writeControl", "patches (corpo carena); writeControl")
    open(os.path.join(nome, "system/controlDict"), "w", newline="\n").write(cd)
    info = dict(geo, X_LE=X_LE, C=C, TAU=TAU, TTE=TTE, R_LE=R_LE, R_TE=R_TE, X_F=X_F, SP_F=SP_F,
                r_disco=rf, Re_L=Re_L, L_V=L, T0=T0, chiave=CHIAVE)
    json.dump(info, open(os.path.join(nome, "carena.json"), "w"), indent=1)
    rbt = lambda xx: float(np.interp(xx, x, r))
    A = lambda xx: np.pi * (float(carena(xx)[0]) ** 2 - rbt(xx) ** 2)
    print(json.dumps(info, indent=1))
    print(f"aree del condotto: presa {A(geo['xc']):.3e}  ventilatore {A(X_F):.3e}  uscita {A(geo['xte']):.3e}")


if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 7.2e-6)
