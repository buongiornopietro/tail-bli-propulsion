"""Caso OpenFOAM "alla Goldschmied": corpo di rivoluzione con presa (aspirazione
dello strato limite) e scarico (getto assiale lungo la coda) ricavati nella parete.

Bilanci (grandezze cinematiche, rho = 1, U = 1, L = 1, cuneo di 5 gradi):
  forza netta sul veicolo  F = F_pareti(corpo+presa+scarico) + sum_{presa,scarico} phi*U_x
  potenza della pompa      P = -sum_{presa,scarico} phi * p0     (p0 = p + |U|^2/2)
Autopropulsione: si regola la portata Q finché F = 0.

uso: python cfd_gold.py <nome> <chiave_profilo> <x_presa_inizio> <x_presa_fine>
                        <Uj_getto> <modello> [Tu]
"""
import json, os, sys, re
import numpy as np
import cfd_caso, cfd_blockmesh_seg

CUNEO = 5.0 / 360.0


def scarico_geom(x, r, x0, larg):
    """Area del tratto di scarico (corpo intero) e componente assiale media della
    normale uscente dal corpo (per il getto assiale)."""
    xx = np.linspace(x0, x0 + larg, 200)
    rr = np.interp(xx, x, r)
    s = np.gradient(rr, xx)
    ds = np.sqrt(1 + s**2)
    A = np.trapezoid(2 * np.pi * rr * ds, xx)
    nx = np.trapezoid(2 * np.pi * rr * (-s), xx) / A      # n_corpo = (-s, 1)/sqrt(1+s^2)
    return A, nx


def main(nome, chiave, xa, xb, Uj, modello="kOmegaSSTLM", Tu=0.001, qmult=1.0):
    prof = json.load(open("out/cfd_profili.json"))[chiave]
    x, r = np.array(prof["x"]), np.array(prof["r"])
    L = x[-1]; x = x / L; r = r / L
    Re_L = prof["ReL"]
    Q_L = qmult * prof["Q_V"] / L**2                         # portata (modello x qmult)
    # larghezza dello scarico per un getto assiale a Uj con la portata di progetto
    rb = np.interp(xb, x, r)
    larg = 0.002
    for _ in range(20):
        A, nx = scarico_geom(x, r, xb, larg)
        larg *= (Q_L / (Uj * nx * A))
    xc = xb + larg
    A, nx = scarico_geom(x, r, xb, larg)
    im = int(np.argmax(r)); xm = float(x[im])
    seg = [(xm, 110, "((0.3 0.4 6) (0.7 0.6 3))", "corpo"),
           (xa, 70, "((0.7 0.6 1) (0.3 0.4 0.15))", "corpo"),
           (xb, 16, "1", "presa"),
           (xc, 12, "1", "scarico"),
           (1.0, 60, "((0.3 0.5 3) (0.7 0.5 0.3))", "corpo")]
    # file standard (fisica, schemi, campi) poi griglia a segmenti
    cfd_caso.caso(nome, x, r, Re_L, 0.7, modello=modello, Tu=Tu)
    bm, h0 = cfd_blockmesh_seg.blockMeshDict(x, r, Re_L, seg, yplus=0.7)
    cfd_caso.scrivi(nome, "system/blockMeshDict", "dictionary", bm)
    # condizioni al contorno di presa e scarico
    Qw = Q_L * CUNEO
    Ujw = Q_L / (nx * A)
    k0 = 1.5 * Tu**2
    bc = {
        "U": (f"presa {{ type flowRateOutletVelocity; volumetricFlowRate {Qw:.6e}; value uniform (0 0 0); }}",
              f"scarico {{ type fixedValue; value uniform ({Ujw:.6f} 0 0); }}"),
        "p": ("presa { type zeroGradient; }", "scarico { type zeroGradient; }"),
        "k": ("presa { type zeroGradient; }", f"scarico {{ type fixedValue; value uniform {k0:.3e}; }}"),
        "omega": ("presa { type zeroGradient; }", f"scarico {{ type fixedValue; value uniform {k0 * Re_L:.4g}; }}"),
        "nut": ("presa { type calculated; value uniform 0; }", "scarico { type calculated; value uniform 0; }"),
        "ReThetat": ("presa { type zeroGradient; }", "scarico { type zeroGradient; }"),
        "gammaInt": ("presa { type zeroGradient; }", "scarico { type fixedValue; value uniform 1; }"),
    }
    for campo, (bp, bs) in bc.items():
        fn = os.path.join(nome, "0", campo)
        if not os.path.exists(fn):
            continue
        t = open(fn).read()
        t = re.sub(r"(boundaryField\s*\{)", r"\1\n        " + bp.replace("\\", "\\\\") + "\n        " + bs.replace("\\", "\\\\"), t, count=1)
        open(fn, "w", newline="\n").write(t)
    # funzioni: forze su tutte le pareti del veicolo, flussi, pressione totale
    cd = open(os.path.join(nome, "system/controlDict")).read()
    cd = re.sub(r"patches(\s+)\(corpo\);(\s+rho)", r"patches\1(corpo presa scarico);\2", cd, count=1)
    extra = """
        ptot
        {
            type            pressure;
            libs            ("libfieldFunctionObjects.so");
            calcTotal       yes;
            calcCoeff       no;
            rho             rhoInf;
            rhoInf          1;
            pRef            0;
            executeControl  timeStep;
            writeControl    writeTime;
        }
        bil_presa
        {
            type            surfaceFieldValue;
            libs            ("libfieldFunctionObjects.so");
            writeControl    timeStep;
            writeInterval   10;
            writeFields     false;
            patch           presa;
            operation       sum;
            weightField     phi;
            fields          (U total(p));
        }
        Q_presa
        {
            type            surfaceFieldValue;
            libs            ("libfieldFunctionObjects.so");
            writeControl    timeStep;
            writeInterval   10;
            writeFields     false;
            patch           presa;
            operation       sum;
            fields          (phi);
        }
        bil_scarico
        {
            type            surfaceFieldValue;
            libs            ("libfieldFunctionObjects.so");
            writeControl    timeStep;
            writeInterval   10;
            writeFields     false;
            patch           scarico;
            operation       sum;
            weightField     phi;
            fields          (U total(p));
        }
        Q_scarico
        {
            type            surfaceFieldValue;
            libs            ("libfieldFunctionObjects.so");
            writeControl    timeStep;
            writeInterval   10;
            writeFields     false;
            patch           scarico;
            operation       sum;
            fields          (phi);
        }
"""
    i = re.search(r"\n\s*yPlus \{", cd).start()
    cd = cd[:i] + "\n" + extra + cd[i:]
    open(os.path.join(nome, "system/controlDict"), "w", newline="\n").write(cd)
    info = dict(xa=xa, xb=xb, xc=xc, larghezza_scarico=larg, A_scarico=A, nx=nx, Q_L=Q_L,
                Q_cuneo=Qw, Uj=Ujw, Re_L=Re_L, h0=h0, L_V=L)
    json.dump(info, open(os.path.join(nome, "gold.json"), "w"), indent=1)
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    a = sys.argv
    main(a[1], a[2], float(a[3]), float(a[4]), float(a[5]),
         a[6] if len(a) > 6 else "kOmegaSSTLM", float(a[7]) if len(a) > 7 else 0.001,
         float(a[8]) if len(a) > 8 else 1.0)
