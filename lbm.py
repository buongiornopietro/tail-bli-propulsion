"""Galleria del vento 2D: Lattice Boltzmann D2Q9 (TRT) con numpy.

Tutte le forme hanno la stessa altezza frontale H: si misura la resistenza
(drag) con il metodo del momentum-exchange sui link di bounce-back.

uso: python lbm.py <forma> <Re> [passi]
"""
import sys, json, time
import numpy as np

PHI = (1 + 5 ** 0.5) / 2

# ---------------------------------------------------------------- forme
def naca_mask(X, Y, x0, y0, L, H):
    """Goccia: profilo simmetrico NACA 4 cifre, lunghezza L, spessore H."""
    t = H / L
    xi = (X - x0) / L
    ok = (xi >= 0) & (xi <= 1)
    xc = np.clip(xi, 0, 1)
    yt = 5 * t * L * (0.2969 * np.sqrt(xc) - 0.1260 * xc - 0.3516 * xc**2
                      + 0.2843 * xc**3 - 0.1036 * xc**4)  # bordo d'uscita chiuso
    return ok & (np.abs(Y - y0) <= yt)


def ellipse_mask(X, Y, x0, y0, L, H):
    return ((X - x0 - L / 2) / (L / 2)) ** 2 + ((Y - y0) / (H / 2)) ** 2 <= 1


def spiral_mask(X, Y, x0, y0, H, turns=1.0, rot=0.0):
    """Spirale aurea 'piena' (guscio di nautilus): r = a*phi^(2θ/π).

    Regione racchiusa dall'ultimo giro della spirale logaritmica aurea;
    scalata perché l'altezza frontale sia H."""
    b = np.log(PHI) / (np.pi / 2)
    th = np.linspace(-6 * np.pi, 2 * np.pi * turns, 4000)
    r = np.exp(b * (th - 2 * np.pi * turns))
    # poligono: ultimo giro esterno + ritorno lungo il giro interno
    sel = th >= 2 * np.pi * turns - 2 * np.pi
    px, py = r[sel] * np.cos(th[sel] + rot), r[sel] * np.sin(th[sel] + rot)
    # scala all'altezza frontale H e piazza con il bordo sinistro in x0
    s = H / (py.max() - py.min())
    px, py = px * s, py * s
    px += x0 - px.min()
    py += y0 - (py.max() + py.min()) / 2
    from matplotlib.path import Path
    pts = np.column_stack([X.ravel(), Y.ravel()])
    m = Path(np.column_stack([px, py])).contains_points(pts).reshape(X.shape)
    return m, px.max() - px.min()


def make_shape(name, X, Y, x0, y0, H):
    if name == "cerchio":
        return ellipse_mask(X, Y, x0, y0, H, H), H
    if name == "ellisse_aurea":
        return ellipse_mask(X, Y, x0, y0, PHI * H, H), PHI * H
    if name == "spirale_aurea":
        return spiral_mask(X, Y, x0, y0, H)
    if name == "spirale_aurea_180":
        return spiral_mask(X, Y, x0, y0, H, rot=np.pi)
    if name.startswith("goccia_"):
        fr = PHI if name == "goccia_aurea" else float(name.split("_")[1])
        return naca_mask(X, Y, x0, y0, fr * H, H), fr * H
    if name == "rettangolo_aureo":
        return ((X >= x0) & (X <= x0 + PHI * H) & (np.abs(Y - y0) <= H / 2)), PHI * H
    raise ValueError(name)


# ---------------------------------------------------------------- LBM
c = np.array([[0, 0], [1, 0], [0, 1], [-1, 0], [0, -1],
              [1, 1], [-1, 1], [-1, -1], [1, -1]])
w = np.array([4 / 9] + [1 / 9] * 4 + [1 / 36] * 4)
opp = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])


def feq(rho, ux, uy):
    cu = 3 * (c[:, 0, None, None] * ux + c[:, 1, None, None] * uy)
    usq = 1.5 * (ux**2 + uy**2)
    return rho * w[:, None, None] * (1 + cu + 0.5 * cu**2 - usq)


def simulate(solid, nu, U, steps, avg_from=0.5, log=None, sponge=0.0,
             outlet="copy"):
    """Simula il flusso attorno alla maschera `solid`; restituisce
    (Fx, Fy, f) con le forze registrate dopo avg_from*steps.

    sponge: frazione finale del dominio in cui la viscosità cresce fino a
    tau = 1 per smorzare i vortici prima dell'uscita (evita riflessioni)."""
    Ny, Nx = solid.shape
    Y, X = np.mgrid[0:Ny, 0:Nx].astype(float)
    tau = np.full((1, Nx), 3 * nu + 0.5)
    if sponge > 0:
        xs = int(Nx * (1 - sponge))
        ramp = ((np.arange(Nx) - xs).clip(0) / (Nx - xs)) ** 2
        tau = tau + (1.0 - tau) * ramp[None, :]
    wp = (1 / tau).astype(np.float32)
    wm = (1 / (0.5 + 0.25 / (tau - 0.5))).astype(np.float32)   # TRT, Λ = 1/4
    f = feq(np.ones((Ny, Nx)), np.full((Ny, Nx), U),
            1e-2 * U * np.sin(2 * np.pi * X / Nx)).astype(np.float32)
    f_in = feq(np.ones((Ny, 1)), np.full((Ny, 1), U),
               np.zeros((Ny, 1))).astype(np.float32)
    fpost = np.empty_like(f)
    wf = w.astype(np.float32)
    sidx = np.flatnonzero(solid)

    bb = []   # (i, indici piatti dei nodi fluidi con vicino a monte solido)
    for i in range(1, 9):
        src_solid = np.roll(solid, (c[i, 1], c[i, 0]), axis=(0, 1))
        bb.append((i, np.flatnonzero(src_solid & ~solid)))

    pairs = [(1, 3), (2, 4), (5, 7), (6, 8)]
    Fx, Fy = [], []
    t0 = time.time()
    for n in range(steps):
        rho = f.sum(0)
        ux = (f[1] + f[5] + f[8] - f[3] - f[6] - f[7]) / rho
        uy = (f[2] + f[5] + f[6] - f[4] - f[7] - f[8]) / rho
        usq = 1.5 * (ux * ux + uy * uy)
        fpost[0] = f[0] - wp * (f[0] - wf[0] * rho * (1 - usq))
        for i, o in pairs:
            cu = 3 * (c[i, 0] * ux + c[i, 1] * uy)
            wr = wf[i] * rho
            ep = wr * (1 + 0.5 * cu * cu - usq)
            em = wr * cu
            fp = 0.5 * (f[i] + f[o]); fm = 0.5 * (f[i] - f[o])
            dp = wp * (fp - ep); dm = wm * (fm - em)
            fpost[i] = f[i] - dp - dm
            fpost[o] = f[o] - dp + dm

        for i in range(9):
            f[i] = np.roll(fpost[i], (c[i, 1], c[i, 0]), axis=(0, 1))
        fx = fy = 0.0
        fpf = fpost.reshape(9, -1); ff = f.reshape(9, -1)
        for i, idx in bb:
            j = opp[i]
            v = fpf[j, idx]
            ff[i, idx] = v
            s = 2 * float(v.sum())
            fx += s * c[j, 0]; fy += s * c[j, 1]
        f[:, :, 0:1] = f_in
        if outlet == "pressione":   # rho = 1, velocità dalla colonna vicina
            g = f[:, :, -2]
            r2 = g.sum(0)
            f[:, :, -1] = feq(np.ones((Ny, 1)),
                              ((g[1] + g[5] + g[8] - g[3] - g[6] - g[7]) / r2)[:, None],
                              ((g[2] + g[5] + g[6] - g[4] - g[7] - g[8]) / r2)[:, None])[:, :, 0]
        else:
            f[:, :, -1] = f[:, :, -2]
        ff[:, sidx] = wf[:, None]

        if not np.isfinite(fx):             # divergenza: interrompi subito
            Fx.append(np.nan); Fy.append(np.nan)
            break
        if n > steps * avg_from:
            Fx.append(fx); Fy.append(fy)
        if log and n % 5000 == 0:
            print(f"{log} step {n}/{steps} Fx={fx:.4f} "
                  f"({time.time() - t0:.0f}s)", flush=True)

    return np.array(Fx), np.array(Fy), f


def run(name, Re, steps=None, H=30, U=0.1, save_field=True):
    Ny = 6 * H
    Nx = 20 * H
    x0, y0 = 4 * H, Ny / 2 + 0.3  # leggero offset: rompe la simmetria perfetta
    Y, X = np.mgrid[0:Ny, 0:Nx].astype(float)
    solid, L = make_shape(name, X + 0.5, Y + 0.5, x0, y0, H)
    area = solid.sum()

    nu = U * H / Re
    if steps is None:
        steps = int(6 * Nx / U)             # ~6 tempi di attraversamento

    t0 = time.time()
    Fx, Fy, f = simulate(solid, nu, U, steps, log=f"{name} Re={Re}")
    q = 0.5 * U**2
    res = dict(forma=name, Re=Re, H=H, L=float(L), LsuH=float(L / H),
               area=int(area),
               Cd=float(Fx.mean() / (q * H)),                  # su altezza frontale
               Cd_area=float(Fx.mean() / (q * np.sqrt(area))),  # a pari area
               Cl=float(Fy.mean() / (q * H)),
               Cl_osc=float(Fy.std() / (q * H)),
               steps=steps, secondi=round(time.time() - t0))
    if save_field:
        rho = f.sum(0)
        ux = (f[1] + f[5] + f[8] - f[3] - f[6] - f[7]) / rho
        uy = (f[2] + f[5] + f[6] - f[4] - f[7] - f[8]) / rho
        vort = (np.roll(uy, -1, 1) - np.roll(uy, 1, 1)
                - np.roll(ux, -1, 0) + np.roll(ux, 1, 0)) / 2
        np.savez_compressed(f"out/{name}_Re{Re:g}.npz", ux=ux, uy=uy,
                            vort=vort, solid=solid, Fx=Fx, Fy=Fy)
    with open(f"out/{name}_Re{Re:g}.json", "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res), flush=True)
    return res


if __name__ == "__main__":
    name, Re = sys.argv[1], float(sys.argv[2])
    steps = int(sys.argv[3]) if len(sys.argv) > 3 else None
    run(name, Re, steps)
