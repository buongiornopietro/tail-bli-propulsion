"""Versione per rivista (punti D1 e D3 della revisione), solo post-elaborazione.

D1  Grandezze navali dei casi RANS in autopropulsione (disco attuatore ideale, rho = U = 1):
      deduzione di spinta          t     = 1 - D_nudo / T
      velocita' totale nel disco   w_tot = 1 - <u>_disco / U   (contiene la velocita' indotta)
      scia effettiva (identita' di spinta, ITTC 7.5-02-03-01.4): il disco in acqua libera alla
        velocita' V_A con la stessa spinta T e la stessa velocita' nel disco <u> (quindi anche
        la stessa potenza T<u>) soddisfa T = 2 rho A <u> (<u> - V_A):
                                   V_A   = <u> - T / (2 rho A <u>),   w_T = 1 - V_A / U
      rendimento di carena         eta_H = (1 - t) / (1 - w_T)
      rendimento in acqua libera   eta_O = V_A / <u>
      rendimento propulsivo        eta_D = eta_H eta_O = D_nudo U / P  (= P_conv / P, elica ideale)
    A = area effettiva Q / u_ventilatore, la stessa dei riferimenti di pari area.
    Scia nominale: corpo NUDO, media di superficie sul disco (r <= 0.03), sui piani
    x = 0.895, 0.900, 0.905 (il volume del disco).
D3  Modello ridotto (bli.py) sullo STESSO corpo dei RANS (ottimo turbolento a ReV 1e6), con la
    presa nella posizione del disco e la stessa portata del disco RANS, contro il risparmio RANS;
    con e senza la penalita' geometrica, e con Q = Q_BL (tutto lo strato limite ingerito).

uso: python rivista_navale.py      -> out/rivista_navale.json
"""
import json
import numpy as np
import misura_carena as mc
import axi, bli

prof = json.load(open("out/cfd_profili.json"))["conv_1e6f"]
xL = np.array(prof["x"]) / prof["x"][-1]; rL = np.array(prof["r"]) / prof["x"][-1]
rb = lambda xx: float(np.interp(xx, xL, rL))
SPICCHI = 72                       # forze per spicchio di 5 gradi -> anello intero

# ---------------- D1
casi = json.load(open("out/misura_carena.json")) + json.load(open("out/misura_v3.json"))
D1 = []
for c in casi:
    if "NON convergente" in c["caso"]:
        continue
    # il campo "P_conv" dei JSON e' D_nudo * U (limite dell'elica ideale), cioe' D_nudo
    T, D0, u = c["T"], c["P_conv"], c["u_disco"]
    A = c["Q"] / c["velocita_media_vent"]
    t, w_tot = 1 - D0 / T, 1 - u
    V_A = u - SPICCHI * T / (2 * A * u)
    w_T = 1 - V_A
    eta_H, eta_O = (1 - t) / (1 - w_T), V_A / u
    D1.append(dict(caso=c["caso"], A=A, t=t, w_tot=w_tot, w_T=w_T, eta_H=eta_H, eta_O=eta_O,
                   eta_D=eta_H * eta_O, P_su_Pconv=c["P"] / D0))
# scia nominale: corpo nudo, sui piani del disco, fino a r 0.03
dati = mc.carica("cfd/conv_ref")
u_piani = []
for xp in (0.895, 0.900, 0.905):
    rr, ux, _ = mc.sezione(dati, xp, rb(xp), 0.03)
    u_piani.append(np.trapezoid(2 * np.pi * rr * ux, rr) / np.trapezoid(2 * np.pi * rr, rr))
w_nom = 1 - float(np.mean(u_piani))

# ---------------- D3
# il modello ridotto rappresenta il propulsore come presa + scarico a parete (bli.py); vicino alla
# coda (raggio piccolo) presa e scarico larghi non stanno sul corpo e bli.power aggiunge una
# penalita' geometrica (vincolo dell'ottimizzatore, non fisica): la si riporta a parte.
x1, r1 = axi.from_profile(prof["x"], prof["r"])
s = x1[-1]                                   # lunghezza del corpo in unita' di V^(1/3)
Q_L = [c for c in casi if c["caso"] == "disco libero (R=0.03)"][0]["Q"]
Q_V = Q_L * s ** 2
ReV = prof["ReV"]
C_conv = bli.conventional(x1, r1, ReV, "forzata")
D3 = []
for xs in (0.85, 0.88, 0.90):
    d = bli.power(x1, r1, ReV, "forzata", xs, Q_V, detail=True)
    dB = bli.power(x1, r1, ReV, "forzata", xs, d["Q_BL"], detail=True)
    D3.append(dict(xs=xs, Q_su_QBL=Q_V / d["Q_BL"], penalita_su_Cconv=100 * d["penalita"] / C_conv,
                   risparmio_senza_penalita=100 * (1 - (d["CP"] - d["penalita"]) / C_conv),
                   risparmio_con_penalita=100 * (1 - d["CP"] / C_conv),
                   risparmio_QBL_senza_penalita=100 * (1 - (dB["CP"] - dB["penalita"]) / C_conv)))

for r in D1:
    print(f"{r['caso']:42s} t {r['t']:+.3f}  w_tot {r['w_tot']:+.3f}  w_T {r['w_T']:+.3f}  "
          f"eta_H {r['eta_H']:.3f}  eta_O {r['eta_O']:.3f}  eta_D {r['eta_D']:.3f}  (P/Pconv {r['P_su_Pconv']:.3f})")
print(f"scia nominale sul disco (corpo nudo, piani 0.895-0.905): w_n = {w_nom:.3f}  "
      f"(piani: {', '.join(f'{1 - v:.3f}' for v in u_piani)})")
for r in D3:
    print(f"D3: presa a x={r['xs']:.2f}, Q del disco RANS (Q/Q_BL {r['Q_su_QBL']:.2f}): risparmio "
          f"{r['risparmio_senza_penalita']:.1f}% senza penalita', {r['risparmio_con_penalita']:.1f}% con "
          f"(penalita' {r['penalita_su_Cconv']:.1f}% di C_conv); con Q = Q_BL {r['risparmio_QBL_senza_penalita']:.1f}%"
          f"  (RANS disco libero: 8.2%)")
json.dump(dict(D1=D1, w_nominale=w_nom, w_nominale_piani=[1 - v for v in u_piani], D3=D3),
          open("out/rivista_navale.json", "w"), indent=1)
