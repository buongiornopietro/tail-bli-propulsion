"""Casi OpenFOAM per la versione 3 della nota (solo generazione dei file, nessun calcolo).

  1) terza griglia, più rada (fattore 1/1.41 per direzione): con la base e la fitta dà
     tre livelli per la stima di Richardson / GCI.
         cfd/v3_conv_rada, cfd/v3_disco_rada, cfd/v3_carena_rada
  2) secondo numero di Reynolds: stesso corpo (conv_1e6f), Re_L 10 volte più alto
     (Re_V circa 1e7); primo strato ricalcolato per y+ < 1.
         cfd/v3_conv_re7, cfd/v3_disco_re7, cfd/v3_carena_re7
     In Allrun_D usare OMEGA=15 (nu_t/nu = 5 come nei casi base, perché nu è 10 volte
     più piccola).
  3) carenature diverse (griglia base, Re base), variando un parametro alla volta
     rispetto al progetto D:
         cfd/v3_carena_E  corda corta e più sottile (C 0.070, TAU 0.04)
         cfd/v3_carena_F  condotto meno convergente (R_TE 0.0275 invece di 0.0249)
         cfd/v3_carena_G  condotto quasi a sezione costante / leggermente divergente
                          (R_TE 0.0290): il ventilatore dovrebbe vedere aria più lenta

uso: python genera_v3.py [griglia|re|carene|tutti]
"""
import json, re, shutil, sys
import numpy as np
import cfd_caso, cfd_blockmesh, cfd_carena

prof = json.load(open("out/cfd_profili.json"))["conv_1e6f"]
X, R = np.array(prof["x"]), np.array(prof["r"]); L = X[-1]; X, R = X / L, R / L
RE_BASE = prof["ReL"]
D_BASE = 6.9574e-06
BASE_CARENA = dict(X_LE=0.855, C=0.095, TAU=0.05, TTE=0.0006, R_LE=0.0293, R_TE=0.0249)

ZONA = """cellZone
        {
            type            cylinder;
            point1          (0.895 0 0);
            point2          (0.905 0 0);
            radius          0.03;
        }"""


def corpo(nome, Re_L, F):
    """Corpo nudo con la topologia di cfd_blockmesh (celle x F in ogni direzione)."""
    shutil.rmtree(nome, ignore_errors=True)
    cfd_caso.caso(nome, X, R, Re_L, 1.0, modello="kOmegaSST")
    bm, h0, _ = cfd_blockmesh.blockMeshDict(X, R, Re_L, 1.0, n_a=round(120 * F), n_b=round(160 * F),
                                            n_scia=round(160 * F), n_bl=round(60 * F), n_est=round(50 * F))
    cfd_caso.scrivi(nome, "system/blockMeshDict", "dictionary", bm)
    return h0


def aggiungi_disco(nome, T0):
    """Disco attuatore libero come in disco_A (x 0.895-0.905, raggio 0.03) + media di U."""
    cfd_caso.scrivi(nome, "constant/fvModels", "dictionary", f"""
    disco
    {{
        type            semiImplicitSource;
        {ZONA}
        volumeMode      absolute;
        sources {{ U {{ explicit ({T0:.4e} 0 0); implicit 0; }} }}
    }}
    """)
    fn = f"{nome}/system/controlDict"
    cd = open(fn).read()
    extra = f"""
        Udisco
        {{
            type            volFieldValue;
            libs            ("libfieldFunctionObjects.so");
            writeControl    timeStep;
            writeInterval   10;
            writeFields     false;
            {ZONA}
            operation       volAverage;
            fields          (U);
        }}
"""
    i = re.search(r"\n\s*yPlus \{", cd).start()
    open(fn, "w", newline="\n").write(cd[:i] + "\n" + extra + cd[i:])


def carena(nome, F=1.0, Re_L=None, T0=7.2e-6, **geom):
    g = dict(BASE_CARENA, **geom)
    for k, v in g.items():
        setattr(cfd_carena, k, v)
    cfd_carena.FATTORE = F
    cfd_carena.RE_L = Re_L
    shutil.rmtree(nome, ignore_errors=True)
    cfd_carena.main(nome, T0)
    for k, v in BASE_CARENA.items():                  # ripristina il progetto D
        setattr(cfd_carena, k, v)
    cfd_carena.FATTORE, cfd_carena.RE_L = 1.0, None


if __name__ == "__main__":
    cosa = sys.argv[1] if len(sys.argv) > 1 else "tutti"
    if cosa in ("griglia", "tutti"):
        F = 1 / 1.41
        corpo("cfd/v3_conv_rada", RE_BASE, F)
        corpo("cfd/v3_disco_rada", RE_BASE, F); aggiungi_disco("cfd/v3_disco_rada", 7.17e-6)
        carena("cfd/v3_carena_rada", F=F, T0=7.5e-6)
    if cosa in ("re", "tutti"):
        Re7 = 10 * RE_BASE
        T7 = D_BASE * 10 ** -0.2                       # stima iniziale: attrito ~ Re^-0.2
        corpo("cfd/v3_conv_re7", Re7, 1.0)
        corpo("cfd/v3_disco_re7", Re7, 1.0); aggiungi_disco("cfd/v3_disco_re7", 1.03 * T7)
        carena("cfd/v3_carena_re7", Re_L=Re7, T0=1.08 * T7)
    if cosa in ("carene", "tutti"):
        carena("cfd/v3_carena_E", C=0.070, TAU=0.04, T0=7.4e-6)
        carena("cfd/v3_carena_F", R_TE=0.0275, T0=7.5e-6)
        carena("cfd/v3_carena_G", R_TE=0.0290, T0=7.5e-6)
    print("fatto:", cosa)
