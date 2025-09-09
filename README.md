# Coordinate Descent Algorithm for Least Absolute Deviations (LAD) Regression

A R implementation of a coordinate descent algorithm for robust Least Absolute Deviations (LAD) regression, providing an efficient and transparent alternative to linear programming-based solutions.

## 📖 Overview

This repository implements a coordinate descent algorithm for solving Least Absolute Deviations (LAD) regression, also known as L₁ regression. Unlike Ordinary Least Squares (OLS) which is sensitive to outliers, LAD regression minimizes the sum of absolute errors, making it robust to anomalous observations.

The algorithm breaks the multidimensional optimization problem into a sequence of one-dimensional subproblems, each solved exactly using medians and weighted medians. This approach guarantees convergence to the global optimum while being simple to implement without specialized optimization software.

## ✨ Key Features

- **Robust Regression**: Handles outliers effectively compared to OLS
- **Coordinate Descent Framework**: Iteratively updates parameters using closed-form solutions
- **Global Convergence**: Guaranteed convergence to optimal solution
- **High-Dimensional Capable**: Works even when p > n (underdetermined systems)
- **No Dependencies**: Pure Python implementation without requiring LP solvers
- **Interpretable**: Easy to understand and modify for educational purposes

## 🧮 Algorithm

The algorithm follows these steps:

1. **Initialize** parameters (zeros or OLS estimates)
2. **Update intercept**: Compute median of partial residuals
3. **Update coefficients**: For each feature, compute weighted median of ratios
4. **Iterate**: Cycle through all parameters until convergence
5. **Convergence check**: Stop when parameter changes fall below tolerance

## 📊 Performance Highlights

- **Robustness**: Matches QuantReg performance on contaminated data
- **Stability**: Converges reliably even in high-dimensional settings (p > n)
- **Efficiency**: Typically converges in under 20 iterations for well-conditioned problems
- **Accuracy**: Achieves predictive performance comparable to established LP solvers

## 🚀 Repository Structure