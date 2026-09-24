import sys, json, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Re = sys.argv[1] if len(sys.argv) > 1 else "200"
ETI = {"cerchio": "Cerchio", "ellisse_aurea": "Ellisse aurea (1:1,618)",
       "spirale_aurea": "Spirale aurea", "spirale_aurea_180": "Spirale aurea (ruotata 180°)",
       "rettangolo_aureo": "Rettangolo aureo", "goccia_aurea": "Goccia aurea (L/H 1,618)",
       "goccia_3": "Goccia L/H 3", "goccia_4.5": "Goccia L/H 4,5"}

res = []
for fn in glob.glob(f"out/*_Re{Re}.json"):
    res.append(json.load(open(fn)))
res.sort(key=lambda r: r["Cd"])

# --- mappe di vorticità
fig, axs = plt.subplots(len(res), 1, figsize=(10, 1.55 * len(res)))
for ax, r in zip(axs, res):
    d = np.load(f"out/{r['forma']}_Re{Re}.npz")
    v = np.ma.masked_array(d["vort"], d["solid"])
    H = r["H"]
    v = v[:, 3 * H:16 * H]
    ax.imshow(v, cmap="RdBu_r", vmin=-0.02, vmax=0.02, origin="lower")
    s = d["solid"][:, 3 * H:16 * H]
    ax.contourf(s, levels=[0.5, 1], colors="k", origin="lower")
    ax.set_title(f"{ETI[r['forma']]}  —  Cd = {r['Cd']:.3f}", fontsize=9, loc="left")
    ax.axis("off")
fig.suptitle(f"Scia (vorticità) dietro ogni forma, flusso da sinistra, Re = {Re}", fontsize=11)
plt.tight_layout()
plt.savefig(f"out/scie_Re{Re}.png", dpi=110)

# --- barre di resistenza
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
for ax, key, tit in [(axs[0], "Cd", "Resistenza a pari ALTEZZA frontale"),
                     (axs[1], "Cd_area", "Resistenza a pari AREA (volume in 2D)")]:
    rr = sorted(res, key=lambda r: r[key])
    col = ["#c0392b" if ("aure" in r["forma"]) else ("#2471a3" if "goccia" in r["forma"] else "#7f8c8d")
           for r in rr]
    col = ["#d4a017" if r["forma"] == "goccia_aurea" else cc for r, cc in zip(rr, col)]
    ax.barh([ETI[r["forma"]] for r in rr], [r[key] for r in rr], color=col)
    for i, r in enumerate(rr):
        ax.text(r[key], i, f" {r[key]:.3f}", va="center", fontsize=8)
    ax.invert_yaxis()
    ax.set_title(tit, fontsize=10)
    ax.set_xlabel("coefficiente di resistenza (più basso = più veloce)")
    ax.spines[["top", "right"]].set_visible(False)
fig.suptitle(f"Simulazione Lattice Boltzmann 2D, Re = {Re}", fontsize=11)
plt.tight_layout()
plt.savefig(f"out/resistenza_Re{Re}.png", dpi=110)

print(f"{'forma':22s} {'L/H':>5s} {'Cd':>7s} {'Cd_area':>8s} {'Cl':>7s} {'Cl_osc':>7s}")
for r in res:
    print(f"{r['forma']:22s} {r['LsuH']:5.2f} {r['Cd']:7.3f} {r['Cd_area']:8.3f} "
          f"{r['Cl']:7.3f} {r['Cl_osc']:7.3f}")
