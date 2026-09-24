"""Casi con griglia infittita (circa 1,41 volte le celle in ogni direzione) per la
verifica di indipendenza dalla griglia: corpo nudo, disco libero, presa intubata."""
import json, os, re, shutil
import numpy as np
import cfd_caso, cfd_blockmesh, cfd_carena

F = 1.41
prof = json.load(open("out/cfd_profili.json"))["conv_1e6f"]
x, r = np.array(prof["x"]), np.array(prof["r"]); L = x[-1]; x, r = x / L, r / L
Re_L = prof["ReL"]

# corpo nudo e disco libero: stessa topologia del caso originale, più celle
for nome in ("cfd/conv_fine", "cfd/disco_fine"):
    shutil.rmtree(nome, ignore_errors=True)
    cfd_caso.caso(nome, x, r, Re_L, 1.0, modello="kOmegaSST")
    bm, h0, _ = cfd_blockmesh.blockMeshDict(x, r, Re_L, 1.0, n_a=round(120 * F), n_b=round(160 * F),
                                            n_scia=round(160 * F), n_bl=round(60 * F), n_est=round(50 * F))
    cfd_caso.scrivi(nome, "system/blockMeshDict", "dictionary", bm)

# disco libero: stesso disco del caso disco_A (x 0.895-0.905, raggio 0.03)
zona = """cellZone
        {
            type            cylinder;
            point1          (0.895 0 0);
            point2          (0.905 0 0);
            radius          0.03;
        }"""
cfd_caso.scrivi("cfd/disco_fine", "constant/fvModels", "dictionary", f"""
    disco
    {{
        type            semiImplicitSource;
        {zona}
        volumeMode      absolute;
        sources {{ U {{ explicit (7.17e-06 0 0); implicit 0; }} }}
    }}
    """)
cd = open("cfd/disco_fine/system/controlDict").read()
extra = f"""
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
"""
i = re.search(r"\n\s*yPlus \{", cd).start()
open("cfd/disco_fine/system/controlDict", "w", newline="\n").write(cd[:i] + "\n" + extra + cd[i:])

# presa intubata (progetto D) con griglia infittita
cfd_carena.FATTORE = F
shutil.rmtree("cfd/carena_fine", ignore_errors=True)
cfd_carena.main("cfd/carena_fine", 7.0e-6)
