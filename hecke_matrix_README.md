# Hecke Operator Matrix Computation

## Overview

The `hecke_matrix.sage` file provides functions to compute matrix representations of Hecke operators acting on spaces of modular forms.

## Main Function

```sage
hecke_operator_matrix(p, N, k=2, cusp=True)
```

Computes the matrix of the p-th Hecke operator (T_p or U_p) acting on cusp forms of level N and weight k.

### Parameters
- `p`: A prime number (the Hecke operator index)
- `N`: The level (positive integer)
- `k`: The weight (default: 2)
- `cusp`: If True, use cusp forms S_k(N); if False, use full modular forms M_k(N) (default: True)

### Returns
A matrix representing the Hecke operator with respect to the standard basis.

## Usage Examples

```sage
# Load the file
load('hecke_matrix.sage')

# Example 1: T_2 on S_2(11) (1-dimensional space)
M = hecke_operator_matrix(2, 11)
# Returns: [-2]

# Example 2: U_2 on S_2(30) (3-dimensional space)
M = hecke_operator_matrix(2, 30)
# Returns: 3x3 matrix

# Example 3: T_3 on S_2(11)
M = hecke_operator_matrix(3, 11)
# Returns: [-1]

# Example 4: Get multiple Hecke operators
matrices = hecke_algebra_generators(11, [2, 3, 5, 7])
# Returns: dictionary {2: T_2, 3: T_3, 5: T_5, 7: T_7}

# Example 5: Full modular forms instead of cusp forms
M = hecke_operator_matrix(2, 11, cusp=False)
# Returns: 2x2 matrix (includes Eisenstein series)
```

## Additional Functions

### `hecke_operator_matrix_on_basis(p, N, basis, k=2)`
Computes the Hecke operator matrix on a custom basis (useful for subspaces).

### `hecke_algebra_generators(N, primes, k=2, cusp=True)`
Computes Hecke operator matrices for multiple primes at once.

## Notes

- When gcd(p, N) = 1, the function computes T_p (the standard Hecke operator)
- When p | N, the function computes U_p (the Atkin-Lehner operator)
- The matrix acts on column vectors of coefficients with respect to the basis S.basis()
- Currently supports trivial character only (can be extended for non-trivial characters)

## Testing

Run the file directly to see examples:
```bash
sage hecke_matrix.sage
```

Or test individual functions:
```bash
sage -c "load('hecke_matrix.sage'); M = hecke_operator_matrix(2, 11); print(M)"
```
