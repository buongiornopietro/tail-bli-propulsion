"""Deriva iterativa della resistenza nelle ultime 400 iterazioni (richiesta del revisore, v3).
Legge postProcessing/forze/*/forceCoeffs.dat dei casi esportati (D = 0.5 * Cd per cuneo).
uso: python deriva_v3.py caso1 caso2 ...
"""
import glob, sys
import numpy as np

for caso in sys.argv[1:]:
    righe = []
    for f in sorted(glob.glob(f"cfd/{caso}/postProcessing/forze/*/forceCoeffs.dat"),
                    key=lambda s: float(s.replace("\\", "/").split("/")[-2])):
        righe += [l.split() for l in open(f) if l.strip() and not l.startswith("#")]
    if len(righe) < 2:
        print(f"{caso:16s} nessun dato"); continue
    it = np.array([float(r[0]) for r in righe]); D = 0.5 * np.array([float(r[2]) for r in righe])
    k = np.searchsorted(it, it[-1] - 400)
    print(f"{caso:16s} it {it[-1]:6.0f}  D {D[-1]:.5e}  ultime {it[-1] - it[k]:.0f} it: {100 * (D[-1] - D[k]) / D[-1]:+.3f}%")
