"""Grafici del caso OpenFOAM assialsimmetrico e confronto con il modello a pannelli.

uso: python plot_cfd.py <cartella_tempo> <chiave_profilo> <titolo>
La cartella tempo deve contenere C, U, p, nut, yPlus in formato ascii.
"""
import re, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import axi


def leggi_lista(testo, vettore):
    m = re.search(r"(\d+)\s*\(\s*", testo)
    n = int(m.group(1))
    corpo = testo[m.end():]
    if vettore:
        v = np.array(re.findall(r"\(([^()]*)\)", corpo[: corpo.find("\n)")])[:n], dtype=object)
        return np.array([list(map(float, s.split())) for s in v])
    return np.array(corpo.split(")")[0].split()[:n], dtype=float)


def campo(fn, vettore=False, patch=None):
    t = open(fn).read()
    if patch is None:
        i = t.index("internalField")
        seg = t[i:t.index("boundaryField")]
    else:
        i = t.index(patch, t.index("boundaryField"))
        seg = t[i:]
        seg = seg[seg.index("value"):]
    if "nonuniform" not in seg.split(";")[0]:
        val = re.search(r"uniform\s+(\(?[^;]*)", seg).group(1)
        return val
    return leggi_lista(seg, vettore)


if __name__ == "__main__":
    d, chiave, titolo = sys.argv[1], sys.argv[2], sys.argv[3]
    C = campo(f"{d}/C", True); U = campo(f"{d}/U", True); p = campo(f"{d}/p")
    nut = campo(f"{d}/nut")
    Cw = campo(f"{d}/C", True, "corpo")
    yp = campo(f"{d}/yPlus", patch="corpo")
    x = C[:, 0]; r = np.hypot(C[:, 1], C[:, 2]); um = np.linalg.norm(U, axis=1)
    # pressione di parete: cella più vicina a ogni faccia del corpo (p ~ costante
    # attraverso lo strato limite; la patch ha zeroGradient e non salva i valori)
    from scipy.spatial import cKDTree
    _, iw = cKDTree(np.c_[x, r]).query(np.c_[Cw[:, 0], np.hypot(Cw[:, 1], Cw[:, 2])])
    pw = p[iw]
    xw = Cw[:, 0]; o = np.argsort(xw)
    prof = json.load(open("out/cfd_profili.json"))[chiave]
    xb, rb = np.array(prof["x"]), np.array(prof["r"]); Lb = xb[-1]
    xc, rc, ue = axi.surface_velocity(xb, rb)

    fig = plt.figure(figsize=(13, 10))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.2, 1.2, 1])
    # 1) velocità attorno al corpo, piano meridiano (specchiato per leggibilità)
    ax = fig.add_subplot(gs[0, :])
    gx, gr = np.meshgrid(np.linspace(-0.3, 2.2, 700), np.linspace(0, 0.45, 140))
    UG = griddata((x, r), um, (gx, gr), method="linear")
    UX = griddata((x, r), U[:, 0], (gx, gr), method="linear")
    UR = griddata((x, r), (U[:, 1] * C[:, 1] + U[:, 2] * C[:, 2]) / np.maximum(r, 1e-12), (gx, gr), method="linear")
    rb_i = np.interp(gx, xb / Lb, rb / Lb, left=0, right=0)
    dentro = gr < rb_i
    for arr in (UG, UX, UR):
        arr[dentro] = np.nan
    cf = ax.contourf(gx, gr, UG, levels=np.linspace(0, 1.25, 26), cmap="RdYlBu_r", extend="both")
    semi = np.c_[np.full(24, -0.29), np.linspace(0.004, 0.44, 24)]
    ax.streamplot(gx, gr, np.nan_to_num(UX), np.nan_to_num(UR), start_points=semi,
                  color="k", linewidth=0.6, arrowsize=0.7, broken_streamlines=False)
    ax.fill_between(xb / Lb, rb / Lb, 0, color="0.2")
    ax.set_aspect("equal"); ax.set_xlim(-0.3, 2.2); ax.set_ylim(0, 0.45)
    ax.set_xlabel("x/L"); ax.set_ylabel("r/L")
    plt.colorbar(cf, ax=ax, label="|U| / U∞", shrink=0.8)
    ax.set_title("Velocità e linee di corrente nel piano meridiano (asse in basso)", loc="left", fontsize=10)
    # 2) zoom sulla coda: strato limite e scia
    ax = fig.add_subplot(gs[1, 0])
    gx2, gr2 = np.meshgrid(np.linspace(0.6, 1.4, 500), np.linspace(0, 0.12, 200))
    UG2 = griddata((x, r), um, (gx2, gr2), method="linear")
    UG2[gr2 < np.interp(gx2, xb / Lb, rb / Lb, left=0, right=0)] = np.nan
    cf = ax.contourf(gx2, gr2, UG2, levels=np.linspace(0, 1.1, 23), cmap="RdYlBu_r")
    ax.fill_between(xb / Lb, rb / Lb, 0, color="0.2")
    ax.set_xlim(0.6, 1.4); ax.set_ylim(0, 0.12)
    ax.set_title("Zoom: strato limite in coda e scia", loc="left", fontsize=10)
    ax.set_xlabel("x/L"); ax.set_ylabel("r/L")
    # 3) viscosità turbolenta (dove c'è turbolenza)
    ax = fig.add_subplot(gs[1, 1])
    NG = griddata((x, r), nut * 4.952e6, (gx2, gr2), method="linear")
    NG[gr2 < np.interp(gx2, xb / Lb, rb / Lb, left=0, right=0)] = np.nan
    cf = ax.contourf(gx2, gr2, np.log10(np.maximum(NG, 1e-2)), levels=20, cmap="magma")
    plt.colorbar(cf, ax=ax, label="log10(ν_t / ν)")
    ax.fill_between(xb / Lb, rb / Lb, 0, color="0.5")
    ax.set_xlim(0.6, 1.4); ax.set_ylim(0, 0.12)
    ax.set_title("Viscosità turbolenta: strato limite e scia", loc="left", fontsize=10)
    ax.set_xlabel("x/L")
    # 4) pressione sul corpo: CFD vs pannelli
    ax = fig.add_subplot(gs[2, 0])
    ax.plot(xw[o], 2 * pw[o], c="#d35400", lw=2, label="CFD RANS (k-ω SST)")
    ax.plot(xc / Lb, 1 - ue**2, "--", c="#2471a3", lw=1.5, label="pannelli (flusso potenziale)")
    ax.invert_yaxis(); ax.set_xlabel("x/L"); ax.set_ylabel("Cp")
    ax.legend(fontsize=8); ax.grid(alpha=.3)
    ax.set_title("Pressione sulla superficie", loc="left", fontsize=10)
    # 5) y+
    ax = fig.add_subplot(gs[2, 1])
    ax.plot(xw[o], yp[o], c="#17a589")
    ax.axhline(1, ls="--", c="k", lw=0.8)
    ax.set_xlabel("x/L"); ax.set_ylabel("y⁺ primo strato"); ax.grid(alpha=.3)
    ax.set_title(f"Risoluzione di parete: y⁺ max = {yp.max():.2f}", loc="left", fontsize=10)
    fig.suptitle(titolo, fontsize=12)
    plt.tight_layout()
    out = f"out/cfd_{chiave}.png"
    plt.savefig(out, dpi=110)
    print(out, " y+ medio %.2f max %.2f" % (yp.mean(), yp.max()))
