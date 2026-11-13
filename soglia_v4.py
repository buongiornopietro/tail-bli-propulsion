# Regia 11/11/2026 (si' di Pietro 10:39 UTC ai punti 1-5): K* col modello corretto della revisione dell'11/11
# (bli.py: massimo lungo la marcia anche per il deficit di Young; axi.py: e^N di Drela-Giles). Due estremi della fascia:
# PHI_MODE "lineare" (presa che prende un campione medio dello strato limite) e "profilo" (strato interno).
# Regime forzata: semi 0, 1, 2 per ogni Re e per ogni estremo; regime naturale (e^N e Michel), seme 0, estremo lineare.
# File out/soglia_Re*_<regime>[_semeN]_v4{lin,pro}[_michel].json: i file precedenti non vengono letti ne' sovrascritti.
# numpy >= 2 (np.trapezoid); 4 processi, circa 13 ore su un PC a 4 core.
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from multiprocessing import Pool
import soglia, bli


def caso(c):
    bli.PHI_MODE = c[-1]           # ogni processo imposta l'estremo prima di ottimizzare
    soglia.KS = [0.0, 0.025, 0.05, 0.075, 0.10, 0.15, 0.20]   # griglia di K piu' fitta (revisione, punto D0)
    return soglia.caso(c[:-1])


casi = [(Re, "forzata", s, 20, f"v4{m[:3]}", m) for m in ("lineare", "profilo") for s in (0, 1, 2) for Re in soglia.RES]
casi += [(Re, "naturale", 0, 20, "v4lin", "lineare") for Re in soglia.RES]
casi += [(Re, "naturale", 0, 20, "v4lin_michel", "michel", "lineare") for Re in soglia.RES]

if __name__ == "__main__":
    with Pool(4) as pool:
        res = list(pool.imap_unordered(caso, casi))
    for c in casi:
        r = caso(c)
        print(c[4], c[1], f"Re {c[0]:.0e}", f"seme {c[2]}", "K* (%)", round(100 * r["K_star"], 2), flush=True)
    print("FATTO", flush=True)
