# Fractal Spacetime Gravity (FSG)

**Author:** Laurent Marie-Louise  
**Date:** December 2025  
**Version:** 4.0 (Major Revision)  
**DOI:** [10.5281/zenodo.17876353](https://doi.org/10.5281/zenodo.17876353)  
**Status:** Preprint (Theoretical Consolidation & Numerical Proof of Principle)

## 📢 What's New in Version 4.0: Rigorous Foundations
This major update elevates FSG from a phenomenological proposal to a consistent **Effective Field Theory (EFT)**. It addresses critical peer-review feedback by introducing:

* **The Fractional Operator:** Replaced the logarithmic ansatz with a rigorous fractional operator $f(X) \sim (L^2/X)^{1/2}$, mathematically proving the emergence of the $k^{-3}$ propagator required for MONDian dynamics.
* **UV Regularization:** Introduced a formal exponential form factor to guarantee vacuum stability and ghost-freedom in the high-energy regime.
* **Analytical Screening:** Provided an explicit derivation of the Vainshtein screening mechanism ($f \sim \sqrt{r}$), ensuring compatibility with Solar System tests (Cassini).
* **Preliminary Boltzmann Solver:** Replaced conceptual sketches with a **Python-based numerical integrator** of the modified Jeans equation, offering a proof-of-principle for the "Fractal Boost" mechanism in the CMB.

---

## Abstract
We propose a modified theory of gravity based on a Fractal and Non-Local Spacetime Geometry (FSG), obtained from an effective action containing the **fractional** non-local operator **$(L^2/\Box^{-1}R)^{1/2}$**. The effective action induces infrared corrections to General Relativity driven by a dimensional flow.

We show that this specific fractional structure implies a reduction of effective spectral dimension towards $\mathbf{d_S \simeq 2}$ in the IR. This geometric reduction is the rigorous cause of the transition to a modified propagator $\mathbf{G(k) \sim k^{-3}}$, which is the requisite condition for MOND dynamics. This leads to:
1.  Flat galactic rotation curves without dark matter.
2.  An exact **Baryonic Tully-Fisher relation** ($V^4 = G M a_0$), derived analytically.
3.  A natural emergence of the acceleration scale **$a_0 \approx c^2\sqrt{\Lambda/3}$**.

To ensure consistency with Solar System tests, we explicitly derive a **source-dependent screening mechanism** restoring General Relativity at small scales. Cosmologically, preliminary numerical integration confirms that the "Fractal Boost" can sustain acoustic oscillations during recombination.

### 🌀 Galaxy Rotation: Newton vs FSG
Visual demonstration of the theory:
* **Left (Newton):** Without Dark Matter, the rotation velocity drops at the edge (Keplerian decline).
* **Right (FSG):** The fractal propagator maintains a high velocity at the edge, naturally producing **flat rotation curves**.

![Galaxy Rotation Simulation](galaxy_rotation.gif)
## 🚀 Unified Physics Simulation

This repository contains the numerical laboratory validating the theory, located in `simulations/`.

### 1. Galactic Dynamics & Structure Formation
* **File:** `simulations/sim_fsg_complet.py`
* **Physics:**
    1.  **3D FFT Convolution:** Solves the modified Poisson equation on a $64^3$ grid using the derived fractal propagator.
    2.  **Linear Growth Solver:** Solves the ODE for cosmic structure formation.
* **Result:** Proves that FSG naturally generates flat rotation curves and explains the early massive galaxies observed by **JWST** ($z \approx 17$).

### 2. CMB "Acid Test" (Boltzmann Check)
* **File:** `simulations/sim_boltzmann_final.py` (New in v4.0)
* **Physics:** Numerical integration of the linear perturbation equations (Modified Jeans Equation) for the baryon-photon fluid.
* **Result:** Demonstrates that the infrared modification of gravity ($G_{\text{eff}}$) successfully re-amplifies the acoustic oscillations, mimicking Dark Matter potential wells without non-baryonic mass.

---

## Key Results Simulated
* **Rotation Curves:** Flatness emerges from $k^{-3}$ scaling (validated via FFT).
* **JWST Observations:** Early collapse at $z \sim 15-20$ (validated via ODE).
* **Solar System:** $\gamma_{PPN} \to 1$ via screening (validated analytically).
* **CMB:** Acoustic peaks restored via Fractal Boost (validated via numerical integration).

## Installation

You need Python installed with scientific libraries:

```bash
pip install numpy matplotlib scipy
```

## 💻 Usage

### 1. Run the Unified Physics Engine (V3.1)
This is the main simulation that reproduces the key results of the paper (Rotation Curves + JWST):

```bash
python simulations/sim_fsg_complet.py
```

### 2. Run CMB Analysis
To visualize the solution to the CMB Third Peak problem:

```bash
python simulations/sim_cmb.py
```

### 3. Generate Paper PDF
To compile the LaTeX source of the article:

```bash
# Ensure you have a LaTeX distribution installed (TeX Live / MiKTeX)
pdflatex main.tex
```

---

## 📄 Citation

If you use this code or theory in your research, please cite **Version 3.1**:

```bibtex
@misc{marielouise2025fsg,
  author       = {Marie-Louise, Laurent},
  title        = {FRACTAL SPATIO-TEMPORAL GRAVITY (FSG): Theoretical Consolidation and Numerical Proof},
  year         = 2025,
  publisher    = {Zenodo},
  version      = {4.0},
  doi          = {10.5281/zenodo.17876353},
  url          = {[https://doi.org/10.5281/zenodo.17876353](https://doi.org/10.5281/zenodo.17876353)}
}
```
