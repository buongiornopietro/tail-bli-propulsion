# Changes to the technical note

## Code release v1.3.0 (journal version, 2026-10-10, doi:10.5281/zenodo.23281153)

New script `rivista_navale.py` (output `out/rivista_navale.json`) for the version submitted to Ocean
Engineering: thrust deduction, effective wake fraction by thrust identity for an ideal actuator disc,
hull, open-water and propulsive efficiencies of the self-propelled RANS cases (free disc: t 0.030,
w_T 0.214, eta_H 1.23, eta_O 0.88, eta_D 1.09); nominal wake over the disc volume (0.24); reduced-order
model on the RANS body with the RANS disc flow rate (9.4-11.1% against 8.2% in RANS) and with the whole
boundary layer ingested (13.8-15.4%). Post-processing only; no RANS result changed.

## Version 3 (2026-10-10, preprint doi:10.5281/zenodo.23272235)

Relative to preprint v2 (doi:10.5281/zenodo.22955086). New RANS cases computed on a second PC with the
same OpenFOAM 14 settings (a v2 case repeated there reproduced the drag to all printed digits).

1. Comparison with the Delft experiment (Della Corte et al. 2021), new Sec. "Comparison with the Delft
   experiment" and table: mechanical energy deposited behind the bare body, same definition as their
   Eq. 10 (kinetic terms plus pressure work, over the momentum-deficit power): 9.9% (RANS, Re_L 4.95e6)
   against about 13% (PIV, Re_L about 2.4e6). Script `delft_confronto.py`.
2. Three more cowls (shorter/thinner; less convergent exit; nearly constant section), new table. The
   shorter, thinner cowl lowers the ducted power from +9.2% to +4.3% (ideal propeller) and from -2.1%
   to -6.2% (equal-area propeller); opening the exit makes it worse. Scripts `genera_v3.py`, `misura_v3.py`.
3. Third (coarser) grid: grid error of the bare-body drag at most 0.5% (GCI, Celik et al. 2008); iterative
   drift of each grid stated (`deriva_v3.py`). The ducted case diverged on the coarse grid and is not reported.
4. Second Reynolds number (Re_L 4.95e7): bare-body drag -29%; free-disc saving 6.7% (ideal propeller) and
   11.9% (equal-area), against 8.2% and 15.3%. The ducted case did not converge at this Reynolds number
   (also with a damped thrust update, `Allrun_smorzato`) and is not reported.
5. Abstract, discussion, conclusions and limitations updated accordingly; consistency fixes from an
   internal review (effective disc areas, ranges, captions, inflow-turbulence exceptions).

No v2 result changed.

## Version 2 (2026-09-25)

Relative to preprint v1 (doi:10.5281/zenodo.22952930).

1. Abstract and conclusions: drag ratios of the golden-ratio shapes now match Table 1
   (2.4–3.6 times in axisymmetric turbulent flow, 1.9–2.5 in the 2D screening; sphere 3.8–7.5).
2. Abstract: saving of the free disc against a propeller of equal disc area stated as 15%
   (v1: "15–18%"; the 17.8% was against a smaller propeller).
3. Table 1, natural-transition rows: the transition criterion is now stated (Michel's criterion,
   Re_x ≤ 5×10^6). Re-evaluated with the e^N criterion, the Re_V = 10^6 body has C_DV = 0.0180
   instead of 0.0072; this is discussed in the text and linked to the nominal e^N shape of the
   robustness section. `axi_opt.py` now sets the criterion per scenario, so it reproduces Table 1.
4. Methods: the boundary-layer closures, transition and separation criteria, shape
   parameterisation and optimiser are described explicitly.
5. References: 13 added (prior axisymmetric shape optimisation and the methods used).

No computed result changed apart from the clarification in item 3.
