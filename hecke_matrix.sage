from sage.all import *

"""
Hecke operator matrix computation for modular forms.

This module provides functions to compute the matrix representation of Hecke operators
acting on spaces of cusp forms (or modular forms) of a given level and weight.
"""

def hecke_operator_matrix(p, N, k=2, cusp=True):
    """
    Compute the matrix of the p-th Hecke operator acting on the space of cusp forms.
    
    This function computes the matrix representation of either T_p (when gcd(p,N)=1)
    or U_p (when p|N) acting on the space S_k(N) of cusp forms of weight k and level N.
    
    Args:
        p: A prime number (the Hecke operator index)
        N: The level (positive integer)
        k: The weight (default: 2)
        cusp: If True, use cusp forms; if False, use full modular forms (default: True)
    
    Returns:
        A matrix over the appropriate base ring representing the Hecke operator.
        The matrix acts on the left on column vectors of coefficients with respect
        to the basis returned by S.basis().
    
    Examples:
        >>> # Hecke T_2 on S_2(11)
        >>> M = hecke_operator_matrix(2, 11)
        >>> M
        [-2]
        
        >>> # Hecke U_2 on S_2(22)
        >>> M = hecke_operator_matrix(2, 22)
        >>> M
        [ 0  1]
        [-1 -1]
        
        >>> # Hecke T_3 on S_2(11)
        >>> M = hecke_operator_matrix(3, 11)
        >>> M
        [-1]
        
        >>> # Full modular forms instead of cusp forms
        >>> M = hecke_operator_matrix(2, 11, k=2, cusp=False)
        >>> M.dimensions()
        (2, 2)
    
    Notes:
        - When gcd(p, N) = 1, this computes the matrix of the Hecke operator T_p
        - When p | N, this computes the matrix of the Hecke operator U_p
        - The matrix is computed with respect to the basis S.basis() where S is the space
        - For trivial character only; for non-trivial characters, extend this function
    """
    if not is_prime(p):
        raise ValueError(f"p = {p} must be a prime number")
    
    if N <= 0:
        raise ValueError(f"Level N = {N} must be a positive integer")
    
    if k <= 0:
        raise ValueError(f"Weight k = {k} must be a positive integer")
    
    # Create the appropriate space
    Space = CuspForms if cusp else ModularForms
    S = Space(N, k)
    
    # Get the Hecke operator matrix
    # Sage automatically uses T_p when gcd(p,N)=1 and U_p when p|N
    T = S.hecke_matrix(p)
    
    return T


def hecke_operator_matrix_on_basis(p, N, basis, k=2):
    """
    Compute the matrix of the p-th Hecke operator on a given basis of forms.
    
    This is useful when you want to compute the Hecke operator on a subspace
    (e.g., a newform space, an eigenspace, or an intersection of spaces).
    
    Args:
        p: A prime number (the Hecke operator index)
        N: The level (positive integer)
        basis: A list of modular forms forming a basis for the subspace
        k: The weight (default: 2)
    
    Returns:
        A matrix representing the Hecke operator on the given basis.
        The (i,j)-th entry is the coefficient of basis[j] in T_p(basis[i]).
    
    Examples:
        >>> S = CuspForms(30, 2)
        >>> basis = S.newforms()[0].basis()  # Basis for a newform space
        >>> M = hecke_operator_matrix_on_basis(7, 30, basis)
        >>> M  # Should be a scalar matrix (eigenvalue)
    
    Notes:
        - All forms in the basis must be in the same space (same level and weight)
        - The matrix is computed by applying the Hecke operator to each basis element
          and expressing the result as a linear combination of the basis
    """
    if not basis:
        raise ValueError("Basis must be non-empty")
    
    if not is_prime(p):
        raise ValueError(f"p = {p} must be a prime number")
    
    # Get the ambient space
    S = basis[0].parent()
    if not hasattr(S, 'level'):
        raise TypeError("Basis elements must be modular forms with a defined level")
    
    # Verify all basis elements are in the same space
    for f in basis[1:]:
        if f.parent() != S:
            raise ValueError("All basis elements must be in the same space")
    
    # Get the Hecke operator
    Tp = S.hecke_operator(p)
    
    # Compute the matrix
    # Apply Tp to each basis element and express as linear combination
    n = len(basis)
    prec = max(int(S.sturm_bound()) + 10, 20)
    
    # Build coefficient matrix for the basis
    basis_coeffs = []
    for f in basis:
        qexp = f.q_expansion(prec)
        basis_coeffs.append([qexp[i] for i in range(prec)])
    
    basis_matrix = matrix(basis_coeffs).transpose()
    
    # Apply Hecke operator and solve
    result_matrix = []
    for f in basis:
        Tpf = Tp(f)
        Tpf_qexp = Tpf.q_expansion(prec)
        Tpf_coeffs = vector([Tpf_qexp[i] for i in range(prec)])
        
        # Solve: basis_matrix * x = Tpf_coeffs
        x = basis_matrix.solve_right(Tpf_coeffs)
        result_matrix.append(x)
    
    return matrix(result_matrix).transpose()


def hecke_algebra_generators(N, primes, k=2, cusp=True):
    """
    Compute matrices for multiple Hecke operators to generate the Hecke algebra.
    
    Args:
        N: The level (positive integer)
        primes: A list of primes (or an integer n to use the first n primes)
        k: The weight (default: 2)
        cusp: If True, use cusp forms; if False, use full modular forms (default: True)
    
    Returns:
        A dictionary mapping each prime p to its Hecke operator matrix.
    
    Examples:
        >>> # Get T_2, T_3, T_5 on S_2(11)
        >>> matrices = hecke_algebra_generators(11, [2, 3, 5])
        >>> matrices[2]
        [-2]
        >>> matrices[3]
        [-1]
        
        >>> # Get the first 10 Hecke operators on S_2(30)
        >>> matrices = hecke_algebra_generators(30, 10)
        >>> len(matrices)
        10
    
    Notes:
        - If primes is an integer n, uses the first n primes
        - Useful for studying the Hecke algebra structure
    """
    # Handle the case where primes is an integer
    if isinstance(primes, (int, Integer)):
        n = primes
        primes = list(primes_first_n(n))
    
    # Compute the matrix for each prime
    result = {}
    for p in primes:
        result[p] = hecke_operator_matrix(p, N, k=k, cusp=cusp)
    
    return result


# Example usage and tests
if __name__ == "__main__":
    print("hi!")
    M = hecke_operator_matrix(3, 22, k=2, cusp=True)
    print(M)
    print(M.charpoly())

print("hoi!")
M = hecke_operator_matrix(11, 60, k=2, cusp=True)
print(M)
print(M.charpoly())
print(M.charpoly().factor())
print(M.jordan_form())

