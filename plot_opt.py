"""Grafici dell'ottimizzazione: convergenza, forma vincente, scia."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import opt
from lbm import simulate

PHI = (1 + 5 ** 0.5) / 2

if __name__ == "__main__":
    d = json.load(open("out/opt_storia.json"))
    hist = [h for h in d["history"] if h["Cd_area"] < 50]
    best = min(hist + ([d["finale"]] if "finale" in d else []),
               key=lambda h: h["Cd_area"])
    p_best = np.array(best["p"])

    fig = plt.figure(figsize=(12, 8.5))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.1])

    # 1) convergenza
    ax = fig.add_subplot(gs[0, 0])
    g = np.array([h["gen"] for h in hist]); cd = np.array([h["Cd_area"] for h in hist])
    ax.scatter(g + np.random.uniform(-.12, .12, len(g)), cd, s=18, c="#95a5a6", label="forme provate")
    gens = sorted(set(g))
    ax.plot(gens, [cd[g == k].min() for k in gens], "o-", c="#2471a3", label="migliore della generazione")
    ax.set_xlabel("generazione"); ax.set_ylabel("Cd a pari area")
    ax.set_title("Convergenza dell'ottimizzazione", fontsize=10, loc="left")
    ax.legend(fontsize=8); ax.spines[["top", "right"]].set_visible(False)

    # 2) snellezza provata vs resistenza
    ax = fig.add_subplot(gs[0, 1])
    fr = np.array([h["LsuH"] for h in hist])
    sc = ax.scatter(fr, cd, c=g, cmap="viridis", s=22)
    ax.axvline(PHI, ls="--", c="#c0392b"); ax.text(PHI, cd.max(), " φ = 1,618", color="#c0392b", fontsize=8, va="top")
    ax.axvline(best["LsuH"], ls="-", c="#2471a3"); ax.text(best["LsuH"], cd.max(), f" ottimo {best['LsuH']:.2f}", color="#2471a3", fontsize=8, va="top")
    ax.set_xlabel("snellezza L/H"); ax.set_ylabel("Cd a pari area")
    ax.set_title("Resistenza in funzione della snellezza", fontsize=10, loc="left")
    plt.colorbar(sc, ax=ax, label="generazione"); ax.spines[["top", "right"]].set_visible(False)

    # 3) profilo ottimo vs partenza vs goccia aurea (stessa area)
    ax = fig.add_subplot(gs[1, :])
    for p, lab, col, lw in [(np.array([np.log(3.0), 0, 0, 0]), "partenza (goccia L/H 3)", "#95a5a6", 1.2),
                            (np.array([np.log(PHI), 0, 0, 0]), "goccia aurea L/H 1,618", "#d4a017", 1.2),
                            (p_best, f"OTTIMO  L/H {best['LsuH']:.2f}", "#2471a3", 2.5)]:
        y, f = opt.profile(p)
        L = np.sqrt(opt.A0 / (2 * np.trapezoid(y, opt.xi)))
        ax.plot(opt.xi * L, y * L, c=col, lw=lw, label=lab)
        ax.plot(opt.xi * L, -y * L, c=col, lw=lw)
    ax.set_aspect("equal"); ax.legend(fontsize=8, loc="upper right")
    ax.set_title("Profili a PARI AREA (flusso da sinistra)", fontsize=10, loc="left")
    ax.spines[["top", "right"]].set_visible(False)

    # 4) scia della forma ottima
    solid, L, H = opt.mask_from_params(p_best)
    Fx, Fy, f = simulate(solid, opt.NU, opt.U, opt.STEPS, avg_from=opt.AVG, outlet="pressione")
    rho = f.sum(0)
    ux = (f[1] + f[5] + f[8] - f[3] - f[6] - f[7]) / rho
    uy = (f[2] + f[5] + f[6] - f[4] - f[7] - f[8]) / rho
    vort = (np.roll(uy, -1, 1) - np.roll(uy, 1, 1) - np.roll(ux, -1, 0) + np.roll(ux, 1, 0)) / 2
    ax = fig.add_subplot(gs[2, :])
    ax.imshow(np.ma.masked_array(vort, solid)[:, 80:], cmap="RdBu_r", vmin=-0.02, vmax=0.02, origin="lower")
    ax.contourf(solid[:, 80:], levels=[0.5, 1], colors="k", origin="lower")
    ax.set_title(f"Scia della forma ottima (vorticità)  —  Cd a pari area = "
                 f"{Fx.mean() / (0.5 * opt.U**2 * np.sqrt(opt.A0)):.3f}", fontsize=10, loc="left")
    ax.axis("off")

    fig.suptitle(f"Ottimizzazione automatica della forma — resistenza minima a pari volume, Re = {opt.Re}",
                 fontsize=12)
    plt.tight_layout()
    plt.savefig("out/ottimizzazione.png", dpi=110)
    np.savez("out/forma_ottima.npz", p=p_best, xi=opt.xi, y=opt.profile(p_best)[0])
    print(json.dumps(best, indent=1))
