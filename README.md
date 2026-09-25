# Optimal streamlined bodies and tail-mounted boundary-layer-ingesting propulsion

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22949388.svg)](https://doi.org/10.5281/zenodo.22949388)

Code and data for the technical note
*"Optimal streamlined bodies and the practical limits of tail-mounted boundary-layer-ingesting
propulsion: a reduced-order and RANS study"* (P. Buongiorno, 2026), source in `nota/nota_tecnica.tex`.

The code was developed with extensive assistance from Claude (Anthropic), an AI assistant.
Comments and variable names are in Italian.

## Requirements
- Python 3.13 with the packages in `requirements.txt` (tested with numpy 2.2.2, scipy 1.15.2,
  matplotlib 3.10.0).
- For the RANS cases: OpenFOAM 14 (openfoam.org), run under WSL/Linux with 4 MPI processes.

All scripts are run from the repository root (they read and write `out/` and `cfd/`).

## What reproduces what

Quick checks (seconds, from the saved data). The plotting scripts write Italian labels; with
the environment variable `FIG_EN=1` they also write the English versions `*_en.png` used in
the note (translation table in `en_fig.py`), e.g. `FIG_EN=1 python plot_axi.py`.

| Item in the note | Command | Input | Output |
|---|---|---|---|
| Fig. 1, 2D drag at Re = 200 | `python plot.py 200` | `out/*_Re200.json/.npz` | `out/resistenza_Re200.png`, `out/scie_Re200.png` |
| Fig. 3, axisymmetric optima | `python plot_axi.py` | `out/axi_risultati.json` | `out/ottimo_3D.png` |
| Fig. 4, K* map | `python plot_mappa_finale.py` | `out/soglia_Re*.json`, `out/v3/` (Michel criterion) | `out/mappa_finale.png` |
| Fig. 5, RANS vs reduced model | `python plot_cfd.py cfd/conv_1e6f/risultati conv_1e6f "<title>"` | `cfd/conv_1e6f/` | `out/cfd_conv_1e6f.png` |
| Fig. 6, tail flow | `python plot_carena.py` | `cfd/disco_A`, `cfd/carena_D` | `out/carena_flusso.png` |
| Slot intake loss K | `python misura_K.py cfd/gold_v2 cfd/K_006 cfd/K_024` | `cfd/gold_v2`, `cfd/K_*` | `out/misura_K.json` |
| Free disc / ducted intake, budget table | `python misura_carena.py` | `cfd/disco_A`, `cfd/carena_D` | `out/misura_carena.json` |
| Grid table | `python misura_fine.py` | `cfd/disco_fine`, `cfd/carena_fine2` | `out/misura_fine.json` (see also `out/griglia.txt`) |
| Reference-propeller table | `python confronto_elica.py` | `out/misura_carena.json`, `out/misura_fine.json` | `out/confronto_elica.json` |
| Fig. 2, "guaranteed" shape | `python plot_robusta.py` | `out/robusta.json` | `out/robusta.png` |

Checked on 2026-09-24: the four measurement scripts reproduce the saved JSON files exactly,
and the figures above are regenerated from the saved data.

Full recomputation (hours to days on a 4-core PC):

| Item | Command | Output |
|---|---|---|
| 2D shape optimisation (not used in the note: coarser lattice) | `python opt.py 8`, then `python plot_opt.py` | `out/opt_storia.json`, `out/ottimizzazione.png` |
| 2D drag of each shape | `python lbm.py <shape> 200` | `out/<shape>_Re200.json` |
| Shape table | `python axi_opt.py` | `out/axi_risultati.json` |
| K* table | `python soglia.py` (all Re); `python soglia.py ripeti 1e7 forzata` (repeatability) | `out/soglia_*.json`, `out/soglia_tutti.json` |
| Goldschmied validation | `python valida_goldschmied.py` | `out/valida_goldschmied.json` |
| "Guaranteed" shape | `python robusta.py` | `out/robusta.json` |

RANS cases (OpenFOAM 14):

| Case | Generator | Run script |
|---|---|---|
| bare body `conv_1e6f`, `conv_ref` | `cfd_caso.py` + `cfd_blockmesh.py` (profiles in `out/cfd_profili.json`) | `Allrun_std` |
| free disc `disco_A` | as above, actuator disc in `constant/fvModels` | `Allrun_std` |
| flush slot `gold_v2`, `K_006`, `K_024` | `cfd_gold.py` + `cfd_blockmesh_seg.py` | `Allrun_gold` |
| ducted intake `carena_D` | `cfd_carena.py` | `Allrun_D` / `Allrun_carena` |
| fine grids `conv_fine`, `disco_fine`, `carena_fine2` | `genera_fini.py` | `Allrun_fini`, `Allrun_fini2` |
| turbulence check along the body | `nut_x.py` | — |

The run-script column was reconstructed from the working folder after the fact; the exact
command sequence of each case has not been re-run from scratch.
`cfd/` contains the exported final fields (`risultati/`) and settings of each case, not the full
time history.
Known pitfalls (spurious laminar SST state, numerical settings) are described in the note.

## Licence
Code: MIT (`LICENSE`). Data, figures and text: CC BY 4.0.
