# Changes to the technical note

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
