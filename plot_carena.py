"""Campo di velocità attorno alla coda: disco libero vs propulsore intubato."""
import en_fig  # noqa: F401  (English figures with FIG_EN=1)
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import misura_carena as mc
import cfd_carena as cc

prof = json.load(open("out/cfd_profili.json"))["conv_1e6f"]
xb = np.array(prof["x"]) / prof["x"][-1]; rbb = np.array(prof["r"]) / prof["x"][-1]
X0, X1, R1 = 0.78, 1.04, 0.065
gx, gr = np.meshgrid(np.linspace(X0, X1, 520), np.linspace(0, R1, 260))

fig, ax = plt.subplots(2, 1, figsize=(12, 6.5), sharex=True)
G = json.load(open("cfd/carena_D/carena.json"))
for k in ("X_LE", "C", "TAU", "TTE", "R_LE", "R_TE"):
    setattr(cc, k, G[k])
casi = [("cfd/disco_A", "Disco attuatore libero (nessuna carenatura)", None),
        ("cfd/carena_D", "Propulsore intubato: carenatura con labbro arrotondato", G)]
for a, (caso, tit, g) in zip(ax, casi):
    x, r, U, p = mc.carica(caso)
    m = (x > X0 - 0.02) & (x < X1 + 0.02) & (r < R1 + 0.01)
    pts = np.c_[x[m], r[m]]
    ux = griddata(pts, U[m, 0], (gx, gr), method="linear")
    ur = griddata(pts, U[m, 1], (gx, gr), method="linear")
    dentro = gr < np.interp(gx, xb, rbb)
    if g is not None:
        ri, ro = cc.carena(gx)
        dentro |= (gx >= g["X_LE"]) & (gx <= g["X_LE"] + g["C"]) & (gr > ri) & (gr < ro)
    print(caso, "frazione NaN", np.isnan(ux).mean())
    ux = np.ma.masked_where(dentro | np.isnan(ux), ux); ur = np.ma.masked_where(dentro | np.isnan(ur), ur)
    cf = a.contourf(gx, gr, ux, levels=np.linspace(0, 1.3, 27), cmap="viridis", extend="both")
    a.streamplot(gx[0], gr[:, 0], np.nan_to_num(ux.filled(0)), np.nan_to_num(ur.filled(0)),
                 density=1.6, color="w", linewidth=0.5, arrowsize=0.6)
    a.fill_between(xb, 0, rbb, color="0.55", zorder=5)
    if g is not None:
        xs = np.linspace(g["X_LE"], g["X_LE"] + g["C"], 300)
        ri, ro = cc.carena(xs)
        a.fill_between(xs, ri, ro, color="0.3", zorder=5)
        xd = [g["X_F"] - g["SP_F"] / 2, g["X_F"] + g["SP_F"] / 2]
        a.fill_between(xd, np.interp(xd, xb, rbb), float(cc.carena(g["X_F"])[0]),
                       color="tab:red", alpha=0.35, zorder=4, lw=0)
        a.text(g["X_F"], 0.0125, "ventilatore", color="tab:red", ha="center", fontsize=8, zorder=6)
    else:
        a.fill_between([0.895, 0.905], np.interp([0.895, 0.905], xb, rbb), 0.03,
                       color="tab:red", alpha=0.35, zorder=4, lw=0)
        a.text(0.9, 0.033, "disco", color="tab:red", ha="center", fontsize=8, zorder=6)
    a.set_title(tit, fontsize=10); a.set_ylabel("r / L"); a.set_ylim(0, R1); a.set_xlim(X0, X1); a.set_aspect("equal")
ax[-1].set_xlabel("x / L")
fig.colorbar(cf, ax=ax, label="velocità assiale u / U", shrink=0.8)
fig.savefig("out/carena_flusso.png", dpi=150, bbox_inches="tight")
print("out/carena_flusso.png")
