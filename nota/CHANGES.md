# Changes to the technical note

## Code release v1.4.0 (corrections to the reduced-order model, 2026-11-13)

Found by two independent reviews of the code; the conclusions of the note are qualitatively unchanged, several
numbers change. The previous files are kept; the new results have the suffix `_v4`.

1. `bli.py`, kinetic-energy defect: the model intake is a wall sink that accelerates the flow just ahead of it, and
   Head's closure lost 3-5% of the defect there, which the kinetic-energy integral equation forbids
   (d/ds(r Ue^3 delta_3) >= 0). The largest value along the forebody is now used.
2. `bli.py`, momentum deficit of the boundary layer that is not ingested (Young's wake drag): same correction.
3. `bli.py`, new switch `PHI_MODE`: "lineare" (default; the intake takes a uniform fraction of the boundary layer,
   as before) or "profilo" (the intake takes only the inner layer; defect fractions from a power-law profile).
4. `axi.py`, e^N transition criterion: now Drela and Giles (1987), dn/ds = dn/dRe_theta(H) (m(H)+1)/2 l(H) / theta.
   Before, the actual dRe_theta/ds was used, which is the same on a flat plate but not on a body whose radius
   grows. Flat-plate transition Re_x at Tu 0.05/0.1/0.3/1%: 2.66e6/1.98e6/1.11e6/4.40e5. This changes the
   natural-transition results (shape robustness, `robusta_v4.py`), not the turbulent ones.
5. New main result: break-even duct loss K* with the propulsor added to the optimal conventional body, only intake
   position and flow rate optimised (`stesso_corpo.py`, exact zero of the lower envelope, no random seeds):
   4.8-2.3% (whole boundary layer ingested) and 9.1-4.7% (inner layer only) for Re_V 3e5-1e8; ideal saving
   16.3-12.7%.
6. The body re-optimised together with the propulsor (`soglia_v4.py`, 7 values of K, three seeds; `fascia_K.py`
   takes the mid-point between the interpolated zero and the envelope of the optimised designs, `inviluppo.py`)
   gives K* 6.4-3.3% and 11.7-6.7% and an ideal saving of 20.0-16.5%, but with a short body whose afterbody behind
   the intake closes at 29-36 degrees; the model counts only friction there (no pressure drag, no separation),
   so these values are an upper estimate. At Re_V 1e8 (whole boundary layer) one seed finds a much shorter body
   exploiting the same assumption; the slender family gives 3.3%.
7. Reduced model on the RANS body (`rivista_navale.py`, `rivista_navale.py profilo`): 7.1-8.1% (whole) and
   13.0-13.1% (inner layer only) against 8.2% in RANS; 13.2-13.5% with Q = Q_BL.
8. Goldschmied validation re-run: 0.01725 against 0.02132 (-19%), unchanged.
9. Ducted case on three grids with self-propulsion balanced to |T-D| < 0.1% (`misura_griglie.py`): +8.7, +8.9,
   +8.9% against the ideal propeller (the values quoted in the note, +9.2% for cowl (d), are balanced to 0.3%).

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
