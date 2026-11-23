# Source bibliography and implementation traceability

References supplied in the original request (R01–R04), two PDFs uploaded to this chat (R05–R06), and one additional peer-reviewed review (R07) are catalogued below.

**Scope:** The application is a *single-material*, locally planar RF and simplified structural **screening** optimizer. The studies on wet radomes, ablation and graded/multilayer walls are included for reference but those submodels are **not implemented**.

## Reference list

### R01 — Nose cone design

**Citation:** Wikipedia contributors (Living reference). *Nose cone design*. Wikipedia.

**Source:** https://en.wikipedia.org/wiki/Nose_cone_design

**Distribution:** External link only. The third-party full text is not redistributed inside the package.

**Relevance:** Profile equations for conic, ogive, elliptical, parabolic, power-series, Haack, and related shapes.

**Application status:** Geometry families implemented with stated approximations; not a computational aerodynamic validation.

### R02 — mmWave Radar Radome Design Guide

**Citation:** C. Kumar, H. U. R. Mohammed, and G. Peake (2021). *mmWave Radar Radome Design Guide*. Texas Instruments application report SWRA705.

**Source:** https://www.ti.com/lit/an/swra705/swra705.pdf

**Direct PDF:** https://www.ti.com/lit/an/swra705/swra705.pdf

**Distribution:** External link only. The third-party full text is not redistributed inside the package.

**Relevance:** Single-wall RF dielectric behavior, thickness resonance, angle sensitivity, and test considerations.

**Application status:** Background support for locally planar dielectric-wall RF model; guide is not implemented verbatim.

### R03 — A Novel Instrument for Real-Time Measurement of Attenuation of Weather Radar Radome Including Its Outer Surface. Part I: The Concept

**Citation:** A. Mancini, J. L. Salazar, R. M. Lebrón, and B. L. Cheong (2018). *A Novel Instrument for Real-Time Measurement of Attenuation of Weather Radar Radome Including Its Outer Surface. Part I: The Concept*. Journal of Atmospheric and Oceanic Technology, 35(5).

**Source:** https://doi.org/10.1175/JTECH-D-17-0083.1

**Direct PDF:** https://arrc.ou.edu/~boonleng/pdf/mancini++18-radome-part-i.pdf

**Distribution:** External link only. The third-party full text is not redistributed inside the package.

**Relevance:** Wet/dirty radome measurement methods and external film/reflection effects.

**Application status:** Literature for future environmental RF model; not implemented.

### R04 — Flight Capabilities of High-Speed-Missile Radome Materials

**Citation:** J. B. Kouroupis (1992). *Flight Capabilities of High-Speed-Missile Radome Materials*. Johns Hopkins APL Technical Digest, 13(3).

**Source:** https://secwww.jhuapl.edu/techdigest/Content/techdigest/pdf/V13-N03/13-03-Kouroupis.pdf

**Direct PDF:** https://secwww.jhuapl.edu/techdigest/Content/techdigest/pdf/V13-N03/13-03-Kouroupis.pdf

**Distribution:** External link only. The third-party full text is not redistributed inside the package.

**Relevance:** Historical high-temperature material tradeoffs, erosion, thermal/mechanical constraints and radar aberration.

**Application status:** Qualitative design context only; no flight qualification, aeroheating, or boresight model implemented.

### R05 — Electromagnetic Performance Analysis of Inhomogeneous Radome Walls Considering Temperature and Ablation

**Citation:** Y. Dan and J. Yang (2023). *Electromagnetic Performance Analysis of Inhomogeneous Radome Walls Considering Temperature and Ablation*. Aerospace, 10, 927.

**Source:** https://doi.org/10.3390/aerospace10110927

**Bundled PDF:** [`references/papers/Dan_Yang_2023_Aerospace_Radome_Temperature_Ablation.pdf`](../references/papers/Dan_Yang_2023_Aerospace_Radome_Temperature_Ablation.pdf)

**Relevance:** Transfer/cascade matrix transmission for wall layers, thermal changes, and ablation of graded materials.

**Application status:** Related mathematical foundation for local slab transmission; monolithic-only model, no ablation or graded wall.

### R06 — A Drop Size Distribution (DSD)-Based Model for Evaluating the Performance of Wet Radomes for Dual-Polarized Radars

**Citation:** J. L. Salazar-Cerreño et al. (2014). *A Drop Size Distribution (DSD)-Based Model for Evaluating the Performance of Wet Radomes for Dual-Polarized Radars*. Journal of Atmospheric and Oceanic Technology, 31, 2409–2430.

**Source:** https://doi.org/10.1175/JTECH-D-13-00208.1

**Bundled PDF:** [`references/papers/Salazar_Cerreno_et_al_2014_Wet_Radomes_DSD.pdf`](../references/papers/Salazar_Cerreno_et_al_2014_Wet_Radomes_DSD.pdf)

**Relevance:** Rain droplet distributions, hydrophobicity, dielectric wet-layer modeling and cross-polarized losses.

**Application status:** Literature for future wet-radome model; precipitation effects not implemented.

### R07 — High temperature ceramic radomes (HTCR) – A review

**Citation:** A. Nag, R. Ramachandra Rao, and P. K. Panda (2021). *High temperature ceramic radomes (HTCR) – A review*. Ceramics International, 47(15), 20793–20806.

**Source:** https://doi.org/10.1016/j.ceramint.2021.04.203

**Distribution:** External link only. The third-party full text is not redistributed inside the package.

**Relevance:** Survey of high-temperature RF-ceramic candidates and manufacturing/property tradeoffs.

**Application status:** Background material screening only; numerical material data in materials.json remain illustrative.

## Traceability by computational subsystem

| Code area | Scientific or software reference | Implementation status |
|---|---|---|
| Parametric nose profiles (`core.py`) | R01 | Approximate parametric construction, not drag/CFD optimization |
| Monolithic TE/TM plane-slab RF (`core.py`) | R02, R05 | Reduced-order transfer-matrix-inspired slab model; no antenna pattern or BSE |
| Thermal boundary effects (`core.py`) | R04, R05 | Prescribed external temperature + 1-D wall conduction, not ablative heating |
| Material candidates (`materials.json`) | R02, R04, R07 | Illustrative records; not verified supplier datasets |
| Wet/iced surface effects | R03, R06 | Future work, not implemented |
| Thickness optimization (`core.py`) | [SciPy SLSQP](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-slsqp.html), [PCHIP](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.PchipInterpolator.html) | Software method references |
| STL/STEP geometry (`exporters.py`) | [Trimesh](https://trimesh.org/), [CadQuery](https://cadquery.readthedocs.io/en/latest/) | Export implementations; CAD validation still required |

## Reproducibility and reuse

- The two bundled papers R05 and R06 originated as PDFs supplied by the user in this conversation. Their hashes are recorded in `catalog.json`. The 2023 Aerospace paper declares an open-access CC BY license.
- The three other third-party PDF URLs (R02–R04) are accessible via direct links and can be retrieved for personal research with the optional `download_external_references.py` script. They are not redistributed in the ZIP due to publisher rights.
- A link to a literature source **does not** mean its complete methodology or measured properties are implemented. No linked paper supplies a validated grade/lot property record for the synthetic default material database.
- Verify current publisher links and material property provenance for any design intended beyond conceptual screening.
