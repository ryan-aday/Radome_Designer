# Radome Studio — monolithic radome concept optimizer

A local Streamlit/Python application that **screens** axial nose geometries and single-material wall-thickness distributions for a civil/industrial RF protective radome. The app produces section-level trends, material comparisons and a **watertight STL / optional solid STEP** geometry in millimeters. The methods are deliberately reduced order and **are not sufficient for flight qualification or manufacturing release**.

## 1. Install and run

Use Python 3.11 or 3.12 with a fresh virtual environment. Recommended on Windows (PowerShell):

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
pip install cadquery  # optional, needed for STEP solid output
streamlit run app.py
```

macOS/Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install cadquery
streamlit run app.py
```

When CadQuery is absent, STL export remains available. If installing CadQuery with `pip` fails on your platform, try its documented conda installation instructions.

Run unit tests:

```bash
pip install pytest
pytest -q
```

### Package layout

```text
app.py                 Streamlit inputs, plots, reports, downloads
core.py                Geometry, RF slab transmission, loads, thermal/stress, optimizer
exporters.py           Normal-offset wall generation, STL and STEP exports
materials.json         Seven illustrative, editable material records
tests/test_core.py     RF/geometry/CAD regression tests
example_design.json    User settings for a demonstration run
example_optimized.*    Example STEP and STL outputs (generated with default inputs)
references/catalog.json                  Source bibliography in machine-readable form
references/README.md                     Full scientific bibliography and source-to-model map
references/papers/*.pdf                  Two user-uploaded original research PDFs
references/download_external_references.py  Optional personal-use PDF fetch script
```

## 2. What is included

- **Nose profiles:** conical, tangent-ogive, secant-ogive, elliptical, parabolic, power-series, von Kármán, LV-Haack, bi-conic.
- **Spherical tip:** optional exact spherical nose cap up to a user-specified joining angle, with C1 cubic-Hermite correction into the selected base profile. It is *not* an exact classical spherically blunted ogive solution for all profiles. Nonmonotonic curves are rejected.
- **Dimensional inputs:** nose length, base diameter, spherical-tip radius and type-specific shape coefficients.
- **RF:** one center frequency, fractional bandwidth, local-incidence envelope (or basic axial-normal proxy), two linear polarizations; transmission based on a single complex-permittivity, homogeneous, locally planar slab (air on both sides). Sampling spans 3 frequencies and the selected angle count.
- **Thermal:** externally prescribed skin-temperature variation from tip to base. 1-D steady wall conduction against internal air convection; internal wall temperature depends on wall thickness and material thermal conductivity.
- **Structural:** magnitude of differential pressure plus axial force, transverse shear and bending moment *resultants*, either constant or arbitrary supported functions of axial position. Thin circular-section tube approximations, a selectable stress thermal-envelope proxy, and reduced tensile/ductile material allowable produce a **screening utilization**. The pointed/transition region and places where radius is <=2× the *maximum permitted* thickness are excluded from stress pass/fail, rather than being falsely analyzed as thin shells.
- **Thickness and material:** 4–12 PCHIP control points. SLSQP multi-start aims to minimize mean one-way insertion loss, mass, and curvature of the thickness profile, subject to thickness and structural limits plus a maximum sampled angle/frequency averaged insertion loss. If infeasible, the optimizer can still return a structural-only candidate flagged as RF-infeasible.
- **Material screening:** 7 illustrative records (PTFE, PEEK, polycarbonate, ABS, fused silica, silicon nitride, alumina). They contain density, dielectric constant, loss tangent, thermal conductivity, CTE, E, Poisson's ratio, representative uniaxial strength, derating and service temperature. Values are NOT manufacturer-qualified data.
- **Alternative geometry:** if no material meets screening, choose **Try alternative nose profiles** for a separate coarse profile ranking; switching an option does **not** overwrite the original study until you select it.
- **Files:** CSV for section stresses/temperatures, JSON for scenario, watertight STL and STEP of revolved shell.

## 3. Load expression reference

Axial coordinate `x` is meters from the nose tip to base; `u = x/L` is the dimensionless axial fraction. Each field returns the *section resultant* in the field's SI units. Examples:

```text
N(x) [newtons]  = 1000*u
V(x) [newtons]  = 150*(1-u)**2
M(x) [N*m]      = 30*u*(1-u)
```

Supported operations: `+ - * / ** %`, comparison operators, `sin`, `cos`, `tan`, `sqrt`, `exp`, `log`, `abs`, `minimum`, `maximum`, `where`, `pi`; no Python function calls, file access, attributes, subscripting or imports. Expressions are parsed with a restricted AST; not free-form Python `eval`.

**Section resultant != line load.** For a distributed normal pressure or force per unit length, the user must first integrate the load into the local internal resultant and apply appropriate signs/boundary conditions. Load models represent a *magnitude-based envelope* and do not solve beam or shell equilibrium.

## 4. Model details and assumptions

### 4.1 Parametric outer geometry

With `x=0` the tip and `x=L` the base, `r(L)=R=D/2`:

- Conic: `r=R*x/L`
- Tangent ogive: `rho=(L²+R²)/(2R)`, `r=sqrt(rho²-(L-x)²)+R-rho`
- Ellipse: `r=R*sqrt(2*x/L-(x/L)²)`
- Parabola: `r=R*(2u-K*u²)/(2-K)`
- Power series: `r=R*u^n`
- Haack: `theta=acos(1-2u)`, `r=R*sqrt((theta-sin(2theta)/2+C*sin(theta)^3)/pi)` with `C=0` for Von Kármán and `C=1/3` for LV-Haack
- Secant ogive: a circle through apex/base with user-selected center axial ratio. The exact permitted ratio depends on L/D.
- Bi-conic: piecewise linear radius vs axial position.

### 4.2 Complex permittivity and transmission

`eps_c = eps_r * (1 - j*tan_delta)`

For incident angle `theta`, let `k_z/k_0 = sqrt(eps_c - sin(theta)^2)` and `delta = k_0*(k_z/k_0)*t`. The TE normalized characteristic *admittance* is `q_TE = sqrt(eps_c - sin(theta)^2)`; the TM admittance is `q_TM=eps_c/sqrt(eps_c - sin(theta)^2)`. The dielectric's 2-port characteristic matrix has `A=D=cos(delta)`, `B=j*sin(delta)/q`, and `C=j*q*sin(delta)`. With air admittance `q0`, `S21=2*q0/(q0*A + q0^2*B + C + q0*D)` and `P_trans=|S21|²`; report `IL_dB=-10*log10(P_trans)`.

An explicit check verifies that a lossless air slab has unit transmission at arbitrary thickness and oblique angle. The **area-weighted RF loss** is a screening surrogate: it does not perform electromagnetic ray tracing, antenna pattern convolution, polarimetric cross coupling, or boresight error analysis. The axial-ray proxy uses only the local normal orientation and an additive scan offset, not a full propagation model.

Temperature dependence of dielectric parameters is represented by editable *linear placeholder coefficients*. Obtain measured εr(f,T), tanδ(f,T) tables for a real material; a broadband high-temperature design cannot be qualified using these default coefficients.

### 4.3 Steady temperature

Given prescribed outer surface temperature `T_s(x)`, inside ambient `T_inf`, and inside convection `h`, heat flux is

`q''(x) = (T_s(x) - T_inf) / (t(x)/k + 1/h)`

and inner wall temperature is `T_in = T_s - q''*t/k`. This assumes the entire outer surface temperature is known; it does **not** solve external convective heating, gas properties, radiative transfer, transient gradients, erosion, ablation, delamination or thermoelastic coupling.

### 4.4 Structural stress *proxies*

For a thin ring at local radius `r`, normal wall thickness `t`, axial resultant `N`, transverse resultant `V`, moment `M`, and pressure difference `Δp`:

`σ_ax ≈ N/(2π r t) + |M|/(π r² t)`

`σ_hoop ≈ |Δp| r / t`, `τ ≈ 2|V|/(π r t)`

Thermal envelope: `σ_T ≈ E α [c_R max(T_mid-T_cabin,0) + c_G |T_out-T_in|] /(1-ν)`.

For ductile bulk materials, stress screen uses a von Mises-like combination. For brittle ceramics, a **maximum principal tensile-stress proxy** is used instead; this still ignores statistical flaw populations, fracture toughness, proof testing, thermal shock and Weibull reliability. Temperature-adjusted strength is linearly derated to an **illustrative** fraction at the illustrative service limit and divided by safety factor. This model is **not valid** near a pointed tip, thick shell, abrupt discontinuity, attachment or shell instability.

### 4.5 CAD construction

The exported wall is a section-of-revolution. Outer profile is sampled axially. Its inner face is approximately the point-by-point offset along the inward meridional normal, remapped to axial coordinates. The inner hollow cavity closes near the nose and remains open at the base. Invalid/self-intersecting offsets raise a message instead of producing silent geometry. The STL triangles form a watertight solid; STEP uses CadQuery to revolve the same closed sampled-section polygon. All exported coordinates are **mm**. The part deliberately lacks a flange and seal.

## 5. Known gaps before engineering use

The following require calibrated extensions:

1. Full-wave electromagnetics (e.g. HFSS/CST/openEMS) or at least validated physical-optics ray tracing with an actual antenna, installed position, phase/scan requirements and boresight-error tolerance.
2. Shell theory / FEM for axisymmetric combined loading, boundary conditions at the mounting flange, thermal strain, out-of-plane deflection, local buckling, eigenmodes and shock response.
3. Manufacturing tolerances, pores, coatings, thickness scatter, moisture and in-service ageing; material tests over frequency, temperature and polarization.
4. Nonstationary heating, through-wall transients, anisotropy, surface radiation, thermal shock, erosion and ablation.
5. Wet/iced radome model: raindrop distribution, surface film and rivulets (the supplied Salazar-Cerreño and Mancini papers are references for later implementation; no precipitation correction is claimed in this release).
6. A true multiobjective Pareto front and a global constrained optimizer with verified convergence.

## 6. Scientific references and origin

The original four user-supplied links, the **two user-uploaded PDFs**, and an additional high-temperature ceramic review are explicitly catalogued below. Every entry is also accessible in the **References & papers** section inside the Streamlit app. The source and implementation-status catalog is maintained in [`references/catalog.json`](references/catalog.json) and detailed in [`references/README.md`](references/README.md).

### Bibliography (original links, uploaded papers, and additional literature)

| ID | Document | Online source | Local PDF | Relation to current app |
|---|---|---|---|---|
| R01 | Nose cone design (Living reference) | [Publisher/source](https://en.wikipedia.org/wiki/Nose_cone_design) | Online only | Geometry families implemented with stated approximations; not a computational aerodynamic validation. |
| R02 | mmWave Radar Radome Design Guide (2021) | [Publisher/source](https://www.ti.com/lit/an/swra705/swra705.pdf) · [PDF](https://www.ti.com/lit/an/swra705/swra705.pdf) | Online only | Background support for locally planar dielectric-wall RF model; guide is not implemented verbatim. |
| R03 | A Novel Instrument for Real-Time Measurement of Attenuation of Weather Radar Radome Including Its Outer Surface. Part I: The Concept (2018) | [Publisher/source](https://doi.org/10.1175/JTECH-D-17-0083.1) · [PDF](https://arrc.ou.edu/~boonleng/pdf/mancini++18-radome-part-i.pdf) | Online only | Literature for future environmental RF model; not implemented. |
| R04 | Flight Capabilities of High-Speed-Missile Radome Materials (1992) | [Publisher/source](https://secwww.jhuapl.edu/techdigest/Content/techdigest/pdf/V13-N03/13-03-Kouroupis.pdf) · [PDF](https://secwww.jhuapl.edu/techdigest/Content/techdigest/pdf/V13-N03/13-03-Kouroupis.pdf) | Online only | Qualitative design context only; no flight qualification, aeroheating, or boresight model implemented. |
| R05 | Electromagnetic Performance Analysis of Inhomogeneous Radome Walls Considering Temperature and Ablation (2023) | [Publisher/source](https://doi.org/10.3390/aerospace10110927) · [PDF](https://www.mdpi.com/2226-4310/10/11/927/pdf) | [Dan_Yang_2023_Aerospace_Radome_Temperature_Ablation.pdf](references/papers/Dan_Yang_2023_Aerospace_Radome_Temperature_Ablation.pdf) | Related mathematical foundation for local slab transmission; monolithic-only model, no ablation or graded wall. |
| R06 | A Drop Size Distribution (DSD)-Based Model for Evaluating the Performance of Wet Radomes for Dual-Polarized Radars (2014) | [Publisher/source](https://doi.org/10.1175/JTECH-D-13-00208.1) | [Salazar_Cerreno_et_al_2014_Wet_Radomes_DSD.pdf](references/papers/Salazar_Cerreno_et_al_2014_Wet_Radomes_DSD.pdf) | Literature for future wet-radome model; precipitation effects not implemented. |
| R07 | High temperature ceramic radomes (HTCR) – A review (2021) | [Publisher/source](https://doi.org/10.1016/j.ceramint.2021.04.203) | Online only | Background material screening only; numerical material data in materials.json remain illustrative. |

### PDF availability and provenance

- **Bundled and downloadable within the app:** R05 Dan & Yang (2023) and R06 Salazar-Cerreño et al. (2014), the **two PDFs uploaded in this conversation**. Both are in `references/papers/` within the ZIP.
- **Link-only PDFs:** R02 Texas Instruments (2021), R03 Mancini et al. (2018), and R04 Kouroupis (1992). Direct original-PDF links are present above, inside the app, and in the machine-readable catalog. They were not repackaged because the externally published copies do not state a general redistribution permission. For your own offline research, you may run `python references/download_external_references.py` to retrieve them into `references/external_downloads/`.
- **Web/article links:** R01 Wikipedia and R07 the journal review. R07's full text may require publisher access.

### Material provenance warning

Values in `materials.json` are **screening placeholders**, not a verified, single-lot or grade-specific data set. The RF and thermal papers support modeling ideas and candidate categories; they do not certify the particular default material rows or their temperature-dependent strength limits. Actual properties and allowables must be validated against the chosen supplier grade, frequency, conditioning, and test temperature before engineering use.
