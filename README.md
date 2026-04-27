# Casio MicroPython Numerical Functions

Numerical methods scripts for Casio calculators running MicroPython.

## Design Constraints

This repository is intentionally procedural and minimal.

The target environment (Casio MicroPython) has limited resources and feature support, so the code avoids:

- object-oriented programming (classes and class hierarchies),
- advanced Python constructs that may be unsupported on-device,
- external dependencies.

## What You Get

- standalone Python scripts for numerical algorithms,
- English comments and user-facing messages,
- implementations designed to be readable and calculator-friendly.

## File Overview

- `aprox.py`: Polynomial approximation using normal equations and Gaussian elimination.
- `crout.py`: Tridiagonal linear system solver (Crout algorithm/Thomas algorithm) with fraction arithmetic.
- `herm_int.py`: Hermite interpolation helper that builds repeated nodes and basis terms step by step.
- `horner.py`: Synthetic polynomial division using Horner's method.
- `lagr_int.py`: Lagrange interpolation polynomial construction with detailed intermediate steps.
- `LLT.py`: Cholesky (`LL^T`) decomposition with optional forward/back substitution.
- `LU.py`: LU factorization with partial pivoting, determinant, and optional matrix inversion.
- `newt_int.py`: Newton interpolation via divided differences.
- `nonlineq.py`: Root-finding methods for nonlinear equations (bisection, regula falsi, secant, Newton).
- `nvil_int.py`: Neville interpolation table for evaluating the interpolating polynomial at a point.
- `poly.py`: Long division of polynomials with fraction coefficients.
- `shawtrb.py`: Tabular algorithm for computing normalized polynomial derivatives at a point.

## Compatibility Goal

The primary goal is reliable execution on constrained calculator hardware, not desktop-style architecture.
