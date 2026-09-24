"""Confronto alla pari: potenza del propulsore in coda (RANS) contro un'elica convenzionale
isolata DELLA STESSA AREA, invece dell'elica ideale di area infinita (P = D U).

Riferimenti (corpo nudo, spinta T = D_nudo, elica in aria libera, teoria del disco):
  ideale : P = D U                                (area infinita, nessuna perdita di getto)
  pari   : P = D U (1 + sqrt(1 + C_T)) / 2,  C_T = T / (½ρU²A)
           con A = area del disco/ventilatore del caso BLI (stessa formula di
           bli.conventional con A_disk).
L'area A viene da Q / <u>_sezione (giro completo), come in misura_carena.py; la spinta è per
cuneo di 5°, quindi T_giro = 72 T. Nessuna simulazione nuova: legge out/misura_carena.json
e out/misura_fine.json.

uso: python confronto_elica.py
"""
import json
import numpy as np

CUNEI = 360 / 5


def fattore_disco(T_cuneo, A):
    CT = CUNEI * T_cuneo / (0.5 * A)
    return CT, (1 + np.sqrt(1 + CT)) / 2


def confronta(r):
    A = r["Q"] / r["velocita_media_vent"]              # area anulare del disco (giro completo)
    D = r["P_conv"]                                     # D_nudo per cuneo (U = 1)
    CT, f = fattore_disco(D, A)
    P_pari = D * f
    return dict(caso=r["caso"], A=A, CT_riferimento=CT, fattore=f, P=r["P"],
                P_ideale=D, P_pari=P_pari,
                variazione_ideale=100 * (r["P"] / D - 1),
                variazione_pari=100 * (r["P"] / P_pari - 1))


if __name__ == "__main__":
    righe = [confronta(r) for r in json.load(open("out/misura_carena.json"))]
    righe += [confronta(r) for r in json.load(open("out/misura_fine.json"))]
    for c in righe:
        print(f"{c['caso']:28s} A={c['A']:.3e}  C_T={c['CT_riferimento']:.3f}  "
              f"P/P_ideale-1 = {c['variazione_ideale']:+5.1f}%   "
              f"P/P_pari-1 = {c['variazione_pari']:+5.1f}%")
    # sensibilità alla misura dell'elica di riferimento (griglia base)
    prof = json.load(open("out/cfd_profili.json"))["conv_1e6f"]
    A_front = np.pi * (max(prof["r"]) / prof["x"][-1]) ** 2
    disco, intub = righe[0], righe[1]
    sens = []
    print("\nsensibilità all'area dell'elica di riferimento (griglia base):")
    for nome, A in [("area ventilatore intubato", intub["A"]), ("area disco libero", disco["A"]),
                    ("2 x area disco libero", 2 * disco["A"]), ("metà sezione frontale", A_front / 2),
                    ("sezione frontale del corpo", A_front), ("infinita (ideale)", np.inf)]:
        f = 1.0 if np.isinf(A) else fattore_disco(disco["P_ideale"], A)[1]
        s = dict(riferimento=nome, A=None if np.isinf(A) else A, fattore=f,
                 disco_libero=100 * (disco["P"] / (disco["P_ideale"] * f) - 1),
                 presa_intubata=100 * (intub["P"] / (intub["P_ideale"] * f) - 1))
        sens.append(s)
        print(f"  {nome:28s} fattore {f:.4f}  disco libero {s['disco_libero']:+5.1f}%  "
              f"presa intubata {s['presa_intubata']:+5.1f}%")
    json.dump(dict(casi=righe, sensibilita=sens, A_frontale=A_front),
              open("out/confronto_elica.json", "w"), indent=1)
