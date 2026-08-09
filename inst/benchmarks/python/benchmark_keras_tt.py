#!/usr/bin/env python3
"""Optional Keras / autodiff TT-P-spline comparator (optional comparator).

This is NOT an R package dependency. Implement outside R CMD check.

Goal: same B-spline bases, TT ranks, λ / penalty as TTPsplines ALS/PIRLS,
but optimize with Adam (or L-BFGS) via autodiff — to contrast conditional
linear / weighted ALS structure vs generic gradient optimization.

Status: stub.
"""

def main():
    print("benchmark_keras_tt.py: stub")
    print("Wire TensorFlow/JAX separately; do not import from R DESCRIPTION.")


if __name__ == "__main__":
    main()
