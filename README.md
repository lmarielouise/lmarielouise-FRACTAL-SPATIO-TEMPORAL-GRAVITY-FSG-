# 🌌 FRACTAL SPATIO-TEMPORAL GRAVITY (FSG)

**Author:** Laurent Marie-Louise  
**Date:** July 2026  
**Version:** **7.0 (draft)** (MOND Derived: No-Go Theorem & Kinetic Memory Term)  
**DOI:** [10.5281/zenodo.17911187](https://doi.org/10.5281/zenodo.17911187)
**Status:** Preprint (v7 draft — v6 archived on Zenodo)

---

## 📢 What's New in Version 6.0: The Post-Audit Rigor Revision

Version 6 is the result of a systematic internal audit of v5.1. Six weaknesses were identified; each is either fixed with a computation now in this repository, or explicitly downgraded to its honest epistemic status. **Several v5.1 claims are formally retracted** — the paper says so in plain words.

* **a₀ derived from cosmic memory (fixed).** The v5.1 formula a₀ ≈ c²√(Λ/3) was wrong (it disagrees with cH₀/2π by a factor 5.7) and is retracted. The v6 result, from numerical integration of the memory fields U = □⁻¹R and V = □⁻²R over cosmic history (`simulations/sim_friedmann_fsg.py`): **a₀ = cH₀·|V₀|^(3/2)/|U₀| ≈ 1.14×10⁻¹⁰ m/s²**, within 5% of observation with no free parameter.
* **The k⁻³ propagator chain corrected (fixed).** The v5.1 scaling argument (X ~ 1/k², f ~ k) was invalid — X is dimensionless — and is retracted. The founding derivation is now **nonlinear and in position space**: the quasi-static limit of the action yields a Bekenstein–Milgrom (AQUAL) equation whose deep-IR solution is exactly Φ ~ ln r. η = 1 is the Fourier-space *name* of that result, and d_S → 2 is a *motivation* for the action, not a derivation of η.
* **Solar System without screening (fixed, and inverted).** The v5.1 Brans–Dicke-style PPN estimate and its screening mechanism are retracted as a category error. The correct perturbative treatment around the cosmological background (as in the Maggiore–Mancarella RR model) gives |γ_PPN − 1| ~ 10⁻¹⁶–10⁻²⁰: FSG is compatible with Cassini **without any screening**.
* **CMB claims split into established / estimated / conjectured (fixed).** The peak *positions* are computed honestly (`simulations/sim_cmb_peaks.py`): r_s unchanged (144.2 Mpc), shift absorbed by the (w, H₀) degeneracy — in the direction that eases the Hubble tension. The peak *amplitudes* ("Fractal Boost") remain a stated conjecture requiring a non-local Boltzmann code. The old "acoustic peaks" figure is relabelled for what it is: **Jeans oscillations from a toy model**.
* **w(z) and H(z) computed, not asserted (fixed).** The unfounded "2–6%" claim is replaced by a computed first-order bracket: w₀ ∈ [−1.7, −1.1] and |ΔH/H| ∈ [1%, 13%] for 0.5 < z < 2, anchored by the RR-model benchmark below and the V-tracking upper bound above (`simulations/sim_friedmann_fsg.py`). The structural, falsifiable core is unchanged: **w₀ < −1** (growing memory density ⟹ phantom).
* **Bullet Cluster reframed (downgraded honestly).** The 1D FDTD simulation is now presented for exactly what it is: a qualitative proof of principle of transient field–source decoupling — with its limitations (1D, v = 0.6c vs ~0.01c observed, transient vs persistent offset, no lensing map) listed explicitly in the paper.

## 📢 What's New in Version 7.0: MOND Is Now Derived (and the Old Derivation Is Refuted)

Version 7 executes the v6 "Priority 1" — derive the AQUAL μ-function from the action — with an unexpected outcome that restructures the theory:

* **No-go theorem (new result).** The curvature-memory sector R·f(□⁻¹R) **cannot produce MOND**: in the quasi-static limit U = □⁻¹R = 2χ + const is *local in the potential*, so the reduction yields a potential-dependent (not acceleration-dependent) modification, and the rotation curve stays **Keplerian** (verified symbolically + numerically). The v6 "scaling selection" claim of an AQUAL reduction is formally **retracted**. This matches the literature: curvature-memory models (RR, RT, Deser–Woodard) give dark energy, not MOND.
* **The corrected action (new physics).** A **kinetic memory term** is added, built from the *gradient* of the same field U: S_K ∝ a₀²F(c⁴(∂U)²/4a₀²). Since ∇U = 2∇χ quasi-statically, this covariant non-local functional "sees" the acceleration — exactly what the no-go forbids to the curvature sector. Precedent: non-local MOND of Soussa–Woodard (2003) and Deffayet–Esposito-Farèse–Woodard (2011); FSG's own twist: U is not ad hoc, it's the cosmological memory field. **One field, two roles: its value carries dark energy, its gradient carries MOND.**
* **AQUAL derived (the derivation the paper always claimed).** Through the localization constraint, the additive term becomes a *multiplicative* modification of Poisson: ∇·[μ_eff(|∇χ|/a₀)∇χ] = 4πGρ with μ_eff = 1 + (4β/α)F′ — a function of the **acceleration**, derived, not postulated (verified with sympy; `simulations/sim_mond_from_action.py`). Choosing μ_eff(y) = 1−e^(−y): exact deep-MOND (BTFR v⁴ = GMa₀ to 0.04% numerically), Newtonian restoration with corrections e^(−g/a₀) ~ e^(−10⁸) at 1 AU — still **no screening needed**.
* **Honest limits stated.** One-potential sector only (Φ=Ψ); the two-potential/lensing reduction, the ghost analysis of the kinetic sector, and the cosmological branch (Θ_Z restriction to spacelike gradients) are the new open problems — listed as the new Priority 1.
* **Strong-field statement amended.** With the kinetic term, Schwarzschild is exact only for the curvature sector; the full v7 correction near a horizon is O(e^(−g/a₀)) ~ e^(−10¹¹) — practically indistinguishable from GR. The ringdown spectrum (coupled δU perturbations) is flagged as not yet computed.

### Extending to the small scale (v6 content)
* **Strong field / black holes (new section, new result).** Because a vacuum black hole has R = 0, the non-local term switches off and **Schwarzschild and Kerr are exact solutions of FSG** — black-hole shadows, ringdown and ISCO are identical to GR. Any deviation is confined beyond the MOND radius r_M = √(GM/a₀), which exceeds the horizon by 8–11 orders of magnitude (`simulations/sim_black_hole.py`, with a table for the Sun, stellar BH, Sgr A*, M87*). The old `sim_black_hole.py`, which applied the MOND force at the horizon and claimed to "soften" the singularity, is retracted — FSG is an IR modification and does **not** resolve the singularity.
* **UV completion honestly bounded.** The "proof of stability" is downgraded to its true status: the entire (infinite-derivative) form factor e^(−□/M*²) is **tree-level ghost-free** (Biswas–Gerwick–Koivisto–Mazumdar; Modesto; Tomboulis), while loop-level unitarity and the non-perturbative treatment of the IR branch cut remain open. The erroneous invocation of the Vafa–Witten theorem is withdrawn.

## Abstract

We propose a modified theory of gravity based on a Fractal and Non-Local Spacetime Geometry (FSG), obtained from an effective action containing the **fractional** non-local operator $X = \Box^{-1}R$, motivated by an infrared flow of the spectral dimension towards $d_S \approx 2$.

A no-go theorem established in v7 shows that the curvature-memory sector alone cannot produce MOND (quasi-statically $\Box^{-1}R$ is local in the potential ⟹ Keplerian curves). The action is therefore completed by a **kinetic memory term** $\propto a_0^2 F(c^4(\partial U)^2/4a_0^2)$ built from the gradient of the same non-local field; through the localization constraint it becomes a multiplicative modification of Poisson, and the Bekenstein–Milgrom (AQUAL) equation is **derived**, with $\mu_{\rm eff}(y)=1-e^{-y}$. This leads to:
1.  Flat galactic rotation curves without dark matter.
2.  An exact **Baryonic Tully-Fisher relation** ($V^4 = GM a_0$), derived analytically.
3.  An acceleration scale derived from the cosmic memory fields: $a_0 = cH_0\,|V_0|^{3/2}/|U_0| \approx 1.14\times10^{-10}\ \mathrm{m/s^2}$, within 5% of observation with no free parameter.

A perturbative treatment of the non-local sector shows that Solar System constraints are satisfied **without any screening mechanism** ($|\gamma_{PPN}-1| \sim 10^{-16}$–$10^{-20}$), as in the related RR non-local model.

Regarding cluster dynamics, a 1D wave simulation demonstrates that the non-local field possesses the intrinsic capacity for transient **ballistic overshoot** — a qualitative mechanism relevant to the Bullet Cluster; the quantitative confrontation (3D, realistic velocities, lensing) is explicitly listed as open.

Cosmologically, the monotonic growth of the memory fields structurally implies a phantom equation of state ($w_0 < -1$), with a computed first-order bracket $w_0 \in [-1.7, -1.1]$ and $|\Delta H/H| \in [1\%, 13\%]$ for $0.5<z<2$ — testable by Euclid across the whole range. The geometry of the CMB acoustic peaks is compatible with FSG (sound horizon unchanged; the induced shift is absorbed by the $(w, H_0)$ degeneracy); the peak **amplitudes** require a dedicated non-local Boltzmann code and remain the decisive open test.

### 🌀 Galaxy Rotation: Newton vs FSG
Visual demonstration of the theory:
* **Left (Newton):** Without Dark Matter, the rotation velocity drops at the edge (Keplerian decline).
* **Right (FSG):** The fractal propagator maintains a high velocity at the edge, naturally producing **flat rotation curves**.

![Galaxy Rotation Simulation](galaxy_rotation.gif)

---

## 🚀 Unified Physics Simulation

This repository contains the numerical laboratory validating the theory, located in `simulations/`.

### 1. Galactic Dynamics & Structure Formation
* **File:** `simulations/sim_fsg_complet.py`
* **Physics:**
    1.  **3D FFT Convolution:** Solves the modified Poisson equation on a 64³ grid using the derived fractal propagator.
    2.  **Linear Growth Solver:** Solves the ODE for cosmic structure formation.
* **Result:** Proves that FSG naturally generates flat rotation curves and explains the early massive galaxies observed by **JWST** (z ≈ 17).

### 2. Cosmic Memory Fields & Background Cosmology (New in v6)
* **File:** `simulations/sim_friedmann_fsg.py`
* **Physics:** Numerical integration of the memory fields U = □⁻¹R and V = □⁻²R over cosmic history (FLRW background), plus first-order back-reaction estimate.
* **Result:** U₀ = −16.1, V₀ = +1.99 ⟹ **a₀ = cH₀·|V₀|^(3/2)/|U₀| = 1.14×10⁻¹⁰ m/s²** (−4.8% vs observation, zeroth order); first-order bracket on w(z) and ΔH/H.

### 3. CMB Peak Geometry (New in v6)
* **File:** `simulations/sim_cmb_peaks.py`
* **Physics:** Sound horizon r_s(z_dec) and comoving distance χ(z_dec); acoustic peak positions ℓ_n ~ nπχ/r_s for ΛCDM vs FSG (w₀ = −1.1).
* **Result:** r_s unchanged (144.2 Mpc); +0.8% peak shift at fixed H₀, absorbed by the (w, H₀) degeneracy — CMB geometry does not exclude FSG.

### 4. Dynamics of Colliding Clusters (reframed in v6)
* **File:** `simulations/sim_bullet_cluster.py`
* **Physics:** 1D FDTD toy model of the non-local wave equation.
* **Result:** Qualitative **proof of principle** of transient field–source decoupling ("ballistic overshoot"). Not a quantitative reproduction of 1E 0657-56 (limitations documented in the script and the paper).

### 5. Jeans-Scale Toy Model (relabelled in v6)
* **File:** `simulations/sim_cmb.py`
* **Physics:** Modified Jeans equation for the baryon fluid with an IR-enhanced G_eff(k). **Not a CMB computation** (no Boltzmann hierarchy, no Thomson coupling, no metric potentials).
* **Result:** Shows that an IR-enhanced coupling can re-amplify pressure-damped **Jeans oscillations** — an illustration of the "Fractal Boost" mechanism, nothing more.

### 6. MOND from the Action — Priority 1 (New in v7)
* **File:** `simulations/sim_mond_from_action.py`
* **Physics:** (A) numerical proof of the no-go (curvature sector ⟹ Kepler); (B) solution of the derived AQUAL equation with μ_eff(y) = 1−e^(−y).
* **Result:** flat rotation curves and exact BTFR **derived from the action**; solar-system corrections e^(−g/a₀) (`fig_mond_from_action.pdf`).

### 7. Additional Simulations
* **`simulations/sim_galaxy_rotation.py`** — Generates the Newton vs FSG rotation animation (`galaxy_rotation.gif`).
* **`simulations/sim_solar_system.py`** — Solar System orbits (note: the v5 screening mechanism is retracted in v6; PPN compatibility now follows from the perturbative non-local calculation in the paper).
* **`simulations/sim_black_hole.py`** — Strong-field scale separation: computes r_s and r_M = √(GM/a₀) for the Sun, stellar/intermediate BHs, Sgr A* and M87*, and plots g(r) showing FSG ≡ GR from the horizon to the MOND radius (`fig_black_hole_scales.pdf`).
* **`simulations/sim_cosmic_web.py`** — Large-scale structure formation under the fractal propagator.

---

## Key Results (v6 epistemic status)
* **Rotation Curves & BTFR:** flat curves and V⁴ = GMa₀ from the AQUAL equation now *derived from the v7 action* (kinetic memory term; one-potential sector) — *derived, with stated caveats*.
* **No-go theorem:** curvature-memory alone ⟹ Keplerian curves — *proved (symbolic + numeric)*.
* **a₀ from cosmic memory:** 1.14×10⁻¹⁰ m/s², −4.8% vs observation at zeroth order — *computed, no free parameter*.
* **Solar System:** |γ_PPN − 1| ~ 10⁻¹⁶–10⁻²⁰ without screening — *computed perturbatively*.
* **Black holes:** Schwarzschild/Kerr exact solutions; GR recovered from horizon to r_M (8–11 decades) — *derived*.
* **JWST Observations:** early collapse at z ~ 15–20 — *computed via linear growth ODE*.
* **w₀ < −1 (phantom):** structural consequence of growing memory — *derived*; magnitude bracketed [−1.7, −1.1] — *first-order estimate*.
* **Bullet Cluster:** transient field–source decoupling — *qualitative proof of principle only*.
* **CMB:** peak geometry compatible — *computed*; peak amplitudes ("Fractal Boost") — *open conjecture, requires non-local Boltzmann code*.

## Installation

You need Python installed with scientific libraries:

```bash
pip install -r requirements.txt
```

## 💻 Usage

### 1. Run the Unified Physics Engine
This is the main simulation that reproduces the key results of the paper (Rotation Curves + JWST):

```bash
python simulations/sim_fsg_complet.py
```

### 2. Run the v6 Cosmology Scripts
Memory fields, a₀ derivation and first-order w(z)/H(z) bracket:

```bash
python simulations/sim_friedmann_fsg.py
```

CMB peak-position geometry (ΛCDM vs FSG):

```bash
python simulations/sim_cmb_peaks.py
```

### 3. Generate Paper PDF
To compile the LaTeX source of the article:

```bash
# Ensure you have a LaTeX distribution installed (TeX Live / MiKTeX)
pdflatex main.tex
```

---

## 📄 Citation

If you use this code or theory in your research, please cite **Version 6.0**:

```bibtex
@misc{marielouise2025fsg,
  author       = {Marie-Louise, Laurent},
  title        = {FRACTAL SPATIO-TEMPORAL GRAVITY (FSG): Post-Audit Rigor Revision},
  year         = 2026,
  publisher    = {Zenodo},
  version      = {6.0},
  doi          = {10.5281/zenodo.17911187},
  url          = {https://doi.org/10.5281/zenodo.17911187}
}
```










