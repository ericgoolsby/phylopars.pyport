# phylopars.pyport

A complete algorithm sandbox for developing Rphylopars 2.0. Originally a port from a Python prototype, it now includes full C++ backends with 10-18x speedups.

## Features

- **ML and REML estimation** of phylogenetic (R) and residual (S) covariance matrices
- **Missing data handling** - arbitrary patterns of missing trait values
- **Within-species variation** - multiple observations per species
- **Evolutionary models** - Brownian Motion (BM) and Pagel's lambda
- **Bastide 2021 algorithm** - O(n) likelihood with analytic gradients
- **EM algorithm** - Expectation-Maximization with closed-form M-step
- **Unified optimizer** - Method sequencing (EM → BFGS → etc.)
- **Multiple backends** - R (debugging), RcppArmadillo (recommended), RcppEigen

## Performance (25 tips, 4 traits)

| Operation | R | Armadillo | Speedup |
|-----------|---|-----------|---------|
| Likelihood (100 evals) | 1.3s | 0.09s | **14x** |
| Gradient (100 evals) | 3.7s | 0.22s | **17x** |
| EM (20 iterations) | 0.55s | 0.03s | **18x** |
| Full Optimizer | 3.5s | 0.37s | **10x** |

## Installation

```r
devtools::install(".")
# Or for development:
devtools::load_all(".")
```

## Quick Start

```r
library(phylopars.pyport)
data(example_tree)
data(example_traits)

# Fast ML estimation with C++ backend (recommended)
result <- optimize_phylopars_cpp(
  example_tree, example_traits,
  method = c("EM", "BFGS"),
  backend = "armadillo"
)
print(result$R)  # Phylogenetic rate matrix
print(result$S)  # Residual covariance matrix

# REML estimation (original approach)
result_reml <- phylopars(example_traits, example_tree)
```

## Key Functions

### High-Level API
| Function | Purpose |
|----------|---------|
| `optimize_phylopars_cpp()` | **Recommended** - Full optimizer with C++ backend |
| `phylopars()` | REML estimation (original interface) |

### ML Likelihood & Gradient (Bastide Algorithm)
| Function | Purpose |
|----------|---------|
| `loglik_bastide_cpp()` | ML log-likelihood |
| `bastide_gradient_cpp()` | Analytic gradients for R, S |

### EM Algorithm
| Function | Purpose |
|----------|---------|
| `em_fit_cpp()` | Full EM with C++ backend |
| `em_fit()` | R implementation |

### Lower-Level Functions
| Function | Purpose |
|----------|---------|
| `restricted_loglik()` | REML likelihood (backend-aware) |
| `loglik_full_matrix()` | O(n³) naive implementation (for validation) |

## Backend Selection

```r
# Recommended - RcppArmadillo (fastest, fully functional)
result <- optimize_phylopars_cpp(tree, data, backend = "armadillo")

# RcppEigen (likelihood/gradient only - EM has known bug)
logL <- loglik_bastide_cpp(R, S, tree, data, backend = "eigen")

# Pure R (slow but useful for debugging)
result <- optimize_phylopars(tree, data)  # Uses R backend
```

## Optimization Methods

```r
# EM only - fast, gets close to optimum
result <- optimize_phylopars_cpp(tree, data, method = "EM")

# EM → BFGS (recommended) - EM initializes, BFGS polishes
result <- optimize_phylopars_cpp(tree, data, method = c("EM", "BFGS"))

# Full chain
result <- optimize_phylopars_cpp(tree, data,
  method = c("EM", "BFGS", "Nelder-Mead"))
```

## Algorithm Overview

1. **Bastide 2021 Tree Traversal**
   - Forward pass (tips → root): Accumulate likelihood contributions
   - Backward pass (root → tips): Compute conditional expectations and gradients
   - O(n) complexity, never constructs full n×n covariance matrix

2. **EM Algorithm**
   - E-step: Conditional expectations E[X_node | data] via tree traversal
   - M-step: Closed-form updates for R, S, μ

3. **Gradient-Based Optimization**
   - Analytic gradients from Bastide backward pass
   - BFGS/L-BFGS-B polish the EM solution

## License

GPL (>= 2)
