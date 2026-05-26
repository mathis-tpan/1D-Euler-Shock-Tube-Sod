# README

# 1D Euler Shock Tube Solver — Sod Problem

## Motivation
This project is a self-initiated, ground-up implementation of a 1D compressible Euler solver in C++, developed during spring 2026 as an introduction to scientific computing and numerical methods in CFD. The Sod shock tube was chosen as the reference test case because it admits an exact analytical solution — making rigorous numerical validation possible — while capturing the three key wave structures of compressible flows: rarefaction fan, contact discontinuity, and shock wave. This solver is the first step of a longer roadmap toward reacting flow simulations, targeting a research internship at the von Karman Institute (VKI) in 2026.

## Physics

The Sod shock tube problem is governed by the 1D Euler equations in conservative form, expressing conservation of mass, momentum, and energy across the domain. The exact analytical solution is derived from two key theoretical tools. First, the Rankine-Hugoniot jump conditions — obtained by integrating the conservation laws across the shock discontinuity — relate the flow states on each side of the shock wave. Second, the Riemann invariant along the rarefaction fan provides a continuous relation between velocity and sound speed in the expansion region. The star state is found using the Brent method (scipy.optimize.brentq) — a robust root-finding algorithm combining bisection and quadratic interpolation.

## Numerical Method

The spatial discretisation uses the **Rusanov (Local Lax-Friedrichs)
scheme**, chosen for its robustness in capturing shock discontinuities.
At each cell interface, the numerical flux is computed as the average
of the left and right physical fluxes, stabilised by a dissipation
term proportional to the maximum local wave speed S_max = max(|u|+c).

Time integration uses an explicit first-order Euler scheme, with the
time step controlled by a CFL condition (CFL = 0.5) applied globally
at each iteration.

The solver is implemented in C++17 following a modular architecture:
- `StateVector` — conservative variable container (struct)
- `Cellule` — mesh cell with derived quantity methods (class)
- `NumericalSchemes` — flux solvers: Rusanov ✓, Roe (TODO), HLLC (TODO)
- `Domaine` — simulation orchestrator: mesh, time loop, I/O (class)

## Results

![Validation](validation.png)

The numerical solution shows good agreement with the exact analytical
solution across all three wave structures — rarefaction fan, contact
discontinuity, and shock wave. The slight smearing of the shock and
contact discontinuity is expected for a first-order Rusanov scheme,
which introduces numerical diffusion near sharp gradients.

**L2 error on density : 0.0069**


## Perspectives

- **MUSCL reconstruction** with minmod slope limiter — 2nd order
  spatial accuracy, reduced numerical diffusion at discontinuities
- **Roe scheme** — approximate Riemann solver with wave decomposition
- **HLLC scheme** — industry standard, better contact discontinuity
  resolution than Rusanov
- **2D extension** — structured mesh, split-dimensional approach
- **Reacting flows** — chemical kinetics and thermodynamics
  (target: VKI internship on reacting flows, 2026)
- **Viscosity** — transition toward full Navier-Stokes,
  eventual integration with OpenFOAM

## References

- Sod, G.A. (1978). *A survey of several finite difference methods
  for systems of nonlinear hyperbolic conservation laws*.
  Journal of Computational Physics, 27, 1-31.
- Toro, E.F. (2009). *Riemann Solvers and Numerical Methods for
  Fluid Dynamics*. Springer.
- Anderson, J.D. (2017). *Fundamentals of Aerodynamics*. McGraw-Hill.
- Rusanov, V.V. (1961). *Calculation of interaction of non-steady
  shock waves with obstacles*. USSR Computational Mathematics.

## Project Structure

├── src/
│   ├── main.cpp
│   ├── StateVector.h
│   ├── Cellule.h / .cpp
│   ├── NumericalSchemes.h / .cpp
│   ├── Domaine.h / .cpp
├── validation.py
├── validation.png
└── README.md
