# Regia 12/11/2026: fascia di K* dai file soglia_*_v4*.json con la regola della revisione Fable xhigh (D, punto a):
# a ogni K il C_P minimo fra i semi e C_conv minimo fra i semi -> K_corda (limite superiore, interpolazione lineare
# come soglia.py); K_inv = zero dell'inviluppo delle rette a progetto fisso di TUTTI i semi (limite inferiore).
# Valore = punto medio; incertezza = max(meta' larghezza, meta' min-max dei K_corda dei singoli semi).
import json, os
import numpy as np
from scipy.optimize import brentq
os.chdir(os.path.dirname(os.path.abspath(__file__)))
ETA = 0.85
RES = [3e5, 1e6, 3e6, 1e7, 3e7, 1e8]


def kcorda(K, dlt):
    if dlt[0] >= 0: return 0.0
    for i in range(len(K) - 1):
        if dlt[i] < 0 <= dlt[i + 1]:
            return K[i] + (K[i + 1] - K[i]) * (-dlt[i]) / (dlt[i + 1] - dlt[i])
    return float("inf")


out = {}
for m in ("lin", "pro"):
    for Re in RES:
        ds = [json.load(open(f"out/soglia_Re{Re:.0e}_forzata{'_seme%d' % s if s else ''}_v4{m}.json")) for s in (0, 1, 2)]
        K = np.array([b["K"] for b in ds[0]["bli"]])
        CP = np.array([[b["CP"] for b in d["bli"]] for d in ds]); Q = np.array([[b["Q"] for b in d["bli"]] for d in ds])
        cc = min(d["C_conv"] for d in ds)
        kc = kcorda(K, CP.min(0) - cc)
        U = lambda k: np.min(CP + (k - K) * Q / ETA) - cc
        ki = brentq(U, 0, 0.6) if U(0) < 0 < U(0.6) else float("nan")
        ks = [kcorda(K, CP[i] - d["C_conv"]) for i, d in enumerate(ds)]
        val = (kc + ki) / 2; inc = max((kc - ki) / 2, (max(ks) - min(ks)) / 2)
        out[f"{m} {Re:.0e}"] = dict(K_corda=100 * kc, K_inv=100 * ki, valore=100 * val, incertezza=100 * inc,
                                    semi=[round(100 * x, 2) for x in ks],
                                    guadagno_K0=round(100 * (CP.min(0)[0] / cc - 1), 1))
        print(f"{m} Re {Re:.0e}: K* = {100 * val:.1f} +- {100 * inc:.1f}  (corda {100 * kc:.2f}, inviluppo {100 * ki:.2f}, "
              f"semi {[round(100 * x, 2) for x in ks]}, guadagno a K=0 {100 * (CP.min(0)[0] / cc - 1):.1f}%)")
for Re in RES:
    for suf in ("naturale_v4lin", "naturale_v4lin_michel"):
        d = json.load(open(f"out/soglia_Re{Re:.0e}_{suf}.json"))
        print(f"{suf} Re {Re:.0e}: K* {100 * d['K_star']:.2f}  x_tr {[round(b['x_trans'], 2) for b in d['bli']]}")
json.dump(out, open("out/fascia_K_v4.json", "w"), indent=1)
