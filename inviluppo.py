"""K*: limite inferiore (inviluppo delle rette a progetto fisso) e superiore (corda, come soglia.py).
A progetto fisso CP è lineare in K: CP_i(K) = CP_i(K_i) + (K - K_i) * Q_i / ETA  (PERDITA = "U", q_din = 1/2).
CP*(K) = min sui progetti  -> concava in K -> la corda fra due K della griglia sta sotto CP*: K_corda >= K*.
Ogni progetto è un limite superiore di CP*: lo zero dell'inviluppo delle rette è K_inv <= K*."""
import json, sys, glob
import numpy as np
from scipy.optimize import brentq
ETA = 0.85
def analizza(fn):
    d = json.load(open(fn)); c = d["C_conv"]; b = d["bli"]
    K = np.array([e["K"] for e in b]); CP = np.array([e["CP"] for e in b]); Q = np.array([e["Q"] for e in b])
    U = lambda k: np.min(CP + (k - K) * Q / ETA) - c
    # controllo: una retta di un altro progetto batte l'ottimo trovato a un K della griglia?
    viol = [(float(K[j]), float(100 * (np.min(CP + (K[j] - K) * Q / ETA) / CP[j] - 1))) for j in range(len(K))]
    kinv = brentq(U, 0, 0.6) if U(0) < 0 < U(0.6) else float("nan")
    return d["Re"], 100 * d["K_star"], 100 * kinv, viol, [round(float(q), 4) for q in Q]
for fn in sys.argv[1:]:
    Re, kc, ki, viol, Q = analizza(fn)
    print(f"{fn.split('/')[-1]:44s} Re {Re:.0e}  K_corda {kc:5.2f}  K_inviluppo {ki:5.2f}  larghezza {kc-ki:4.2f}   "
          f"min retta/ottimo-1 (%): {[round(v[1],2) for v in viol]}  Q: {Q}")
