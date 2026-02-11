from sage.all import *

"""
Modular forms toolkit with robust newform support.

All functions automatically handle both regular modular form elements and newforms.
When a newform causes issues (TypeError: "can't initialize vector from nonzero non-list"),
the functions automatically convert via q-expansion or ambient space as needed.

Key functions:
- HeckeT, HeckeU: Hecke operators (work with newforms)
- alpha, beta: Degeneracy maps (work with newforms)
- LinearDependence, Intersection, ContainedIn: Linear algebra on modular forms (work with newforms)
"""

# 1) Spaces at levels N, Np, Np^2 (weight k, trivial character; set cusp=True for cusp forms)
def spaces_N_p(N, p, k=2, cusp=False):
    Space = CuspForms if cusp else ModularForms
    return Space(N, k), Space(N*p, k), Space(N*p*p, k)


# Internal: ensure we were passed a modular form element, not a raw q-expansion
def _space_of_form(f):
    S = f.parent()
    if not hasattr(S, 'level'):
        raise TypeError("Expected a modular form element (e.g. from M.basis()), not a q-expansion. Construct the form element in a space before applying Hecke/Atkin–Lehner/degeneracy.")
    return S

# Internal: convert newforms to regular modular forms via q-expansion
# This avoids "can't initialize vector from nonzero non-list" errors
def _normalize_form(f, target_space=None):
    """
    Convert a modular form (especially newforms) to a regular vector-based form.
    
    Args:
        f: A modular form element
        target_space: Optional target space (e.g., ambient_module()). If None, uses f.parent()
    
    Returns:
        A regular modular form in target_space that can be used with operators
    """
    S = f.parent()
    if target_space is None:
        target_space = S
    
    # Get sufficient precision for the conversion
    prec = max(int(target_space.sturm_bound()) + 10, 20)
    qexp = f.q_expansion(prec)
    return target_space(qexp)

# Hecke T_p (requires p coprime to level)
def HeckeT(f, p):
    S = _space_of_form(f)
    N = S.level()
    if gcd(p, N) != 1:
        raise ValueError("T_p is for p coprime to the level; use HeckeU when p | level.")
    f_normalized = _normalize_form(f, S)
    return S.hecke_operator(p)(f_normalized)

# Hecke U_p (requires p | level)
def HeckeU(f, p):
    S = _space_of_form(f)
    N = S.level()
    if N % p != 0:
        raise ValueError("U_p requires p | level.")
    f_normalized = _normalize_form(f, S)
    return S.hecke_operator(p)(f_normalized)

# Partial Atkin–Lehner W_d (d must be a Hall divisor: d | N and gcd(d, N/d) = 1)
def AtkinLehner(f, d):
    S = _space_of_form(f)
    N = S.level()
    if N % d != 0 or gcd(d, N//d) != 1:
        raise ValueError("d must be a Hall divisor of the level (d | N and gcd(d, N/d)=1).")
    
    # Note: Sage's implementation of Atkin-Lehner for ModularForms is limited.
    # For composite levels, it only works when d = N (the full Atkin-Lehner involution).
    # Partial Atkin-Lehner operators W_d (d != N) are only implemented for ModularSymbols.
    # For prime power levels, all operators should work.
    
    f_normalized = _normalize_form(f, S)
    try:
        return S.atkin_lehner_operator(d)(f_normalized)
    except NotImplementedError:
        N_factorization = factor(N)
        if d != N and len(N_factorization) > 1:
            raise NotImplementedError(
                f"Sage does not implement partial Atkin-Lehner operators W_d for ModularForms "
                f"at composite levels. Level {N} = {N_factorization}. "
                f"You can use AtkinLehner(f, {N}) for the full involution, "
                f"or work with ModularSymbols instead of ModularForms for partial operators."
            ) from None
        else:
            raise

# Degeneracy maps α, β: level N -> level Np (weight preserved)
# α = d=1, β = d=p
def alpha(f, p):
    S = _space_of_form(f)
    ambient = S.ambient_module()
    f_normalized = _normalize_form(f, ambient)
    return ambient.degeneracy_map(ambient.level()*p, 1)(f_normalized)

def beta(f, p):
    S = _space_of_form(f)
    ambient = S.ambient_module()
    f_normalized = _normalize_form(f, ambient)
    return ambient.degeneracy_map(ambient.level()*p, p)(f_normalized)

# ---- Linear algebra helpers on q-expansions up to a Sturm bound ----

def _sturm_prec_for_forms(forms):
    if not forms:
        return 0
    S = forms[0].parent()
    return int(S.sturm_bound())

def _coeff_vector(f, prec):
    """Vector [a_0, a_1, ..., a_{prec-1}] of q-expansion coefficients."""
    qexp = f.q_expansion(prec)
    K = qexp.parent().base_ring()
    return vector(K, [qexp[n] for n in range(prec)])

def _coeff_matrix_columns(forms, prec):
    """Returns a matrix with columns equal to coefficient vectors of the given forms."""
    if not forms:
        raise ValueError("Need at least one form.")
    vecs = [_coeff_vector(f, prec) for f in forms]
    return matrix(vecs[0].parent().base_ring(), vecs).transpose()

def _check_forms_compatible(forms):
    """Check that all forms have the same level and weight."""
    if not forms:
        return
    first_parent = forms[0].parent()
    for f in forms[1:]:
        f_parent = f.parent()
        if (hasattr(first_parent, 'level') and hasattr(f_parent, 'level') and
            hasattr(first_parent, 'weight') and hasattr(f_parent, 'weight')):
            if (first_parent.level() != f_parent.level() or 
                first_parent.weight() != f_parent.weight()):
                raise ValueError("All forms must have the same level and weight.")
        elif f_parent is not first_parent:
            # Allow some flexibility for newforms vs basis elements from the same space
            pass

# 3) Linear dependence among forms: returns (is_independent, kernel_basis)
# kernel_basis is a list of relations [c1,...,cn] with sum c_i f_i = 0
def LinearDependence(forms, prec=None):
    if not forms:
        return False, []
    
    _check_forms_compatible(forms)
    
    if prec is None:
        prec = _sturm_prec_for_forms(forms)
    
    # Rows = forms, so that left_kernel gives relations among rows
    coeff_vectors = [_coeff_vector(f, prec) for f in forms]
    base_ring = coeff_vectors[0].parent().base_ring()
    M = matrix(base_ring, coeff_vectors)
    ker = M.left_kernel()
    return (ker.dimension() != 0, [v for v in ker.basis()])

# 4) Intersection of subspaces U, V given by their bases (lists of forms in the same space)
# Returns a list of nonzero forms forming a basis for the intersection (up to the chosen precision).
def Intersection(U_basis, V_basis, prec=None):
    if not U_basis or not V_basis:
        return []
    
    all_forms = U_basis + V_basis
    _check_forms_compatible(all_forms)
    
    if prec is None:
        prec = _sturm_prec_for_forms(all_forms)

    U_mat = _coeff_matrix_columns(U_basis, prec)  # prec x rU
    V_mat = _coeff_matrix_columns(V_basis, prec)  # prec x rV
    A = U_mat.augment(-V_mat)                     # prec x (rU + rV)

    ker = A.right_kernel().basis()               # vectors (x, y) with U_mat*x = V_mat*y

    rU = len(U_basis)
    first_parent = all_forms[0].parent()
    intersection_basis = []
    for w in ker:
        x = w[:rU]
        g = sum((x[j] * U_basis[j] for j in range(rU)), first_parent(0))
        if g != 0:
            intersection_basis.append(g)
    return intersection_basis

# 5) Containment test f ∈ span(U_basis)
# Returns (contained, coefficients) where coefficients solve sum c_j U_j = f if contained.
def ContainedIn(f, U_basis, prec=None):
    if not U_basis:
        return False, None
    
    all_forms = [f] + U_basis
    _check_forms_compatible(all_forms)
    
    if prec is None:
        prec = _sturm_prec_for_forms(all_forms)

    try:
        U_mat = _coeff_matrix_columns(U_basis, prec)     # prec x r
        v_f  = _coeff_vector(f, prec)                    # length prec
        coeffs = U_mat.solve_right(v_f)                  # solves U_mat * c = v_f
        return True, coeffs
    except (ValueError, ArithmeticError):
        return False, None


for N in range(1, 401):
  print("trying N", N)
  E2N = EisensteinForms(N, 2);
  dim_eis = len(E2N.eisenstein_series())
  print(f"there are {dim_eis} Eisenstein series at this level")

 p = Primes(modulus=N)[0]
  lin_dep, ker = LinearDependence([alpha(h, p) for h in E2N.eisenstein_series()] + [beta(h, p) for h in E2N.eisenstein_series()])

  if lin_dep:
    print("found a collision!", N, ker)

"""
S15, S30 = CuspForms(15, 2), CuspForms(30, 2)
f = S15.newforms()[0]
g = S30.newforms()[0]

#S11, S22 = CuspForms(11, 2), CuspForms(22, 2)
#f = S11.newforms()[0]

print(f)
print(g)
print("---")
print(HeckeT(f, 2))                        # Hecke T_p works with newforms
print(HeckeU(g, 2))                        # Hecke U_p works with newforms
print(LinearDependence([g, HeckeU(g, 2)])) # Linear algebra works with newforms
print(LinearDependence([f, HeckeT(f, 2)]))
print("---")
# should print False
print(HeckeU(beta(f, 2), 2) == alpha(f, 2))
print(HeckeU(beta(g, 2), 2) == alpha(g, 2))
print("---")
# should print True
print(HeckeU(alpha(f, 2), 2) == beta(f, 2))
print(HeckeU(alpha(g, 2), 2) == beta(g, 2))
print("---")
print(LinearDependence([alpha(g,2), beta(g,2)]))

#h = M11.eisenstein_series()[0]
#[(p, h.coefficient(p)) for p in primes_first_n(20)]

# Example usage
# Spaces
M11, M22, M44 = spaces_N_p(11, p=2, k=2, cusp=False)
S11, S22, S44 = spaces_N_p(11, p=2, k=2, cusp=True)

# Pick a form at level 11
f = M11.basis()[0]
# Hecke T_2 (since 2 ∤ 11)
T2f = HeckeT(f, 2)

# Degeneracy maps α, β: level 11 -> level 22
f_alpha = alpha(f, 2)                     # lives in M22
f_beta  = beta(f, 2)                      # lives in M22

# Hecke U_2 at level 22 (since 2 | 22)
g = f_alpha
U2g = HeckeU(g, 2)

# Linear dependence among some forms at level 22
indep, rels = LinearDependence([f_alpha, f_beta, M11.degeneracy_map(22,1)(T2f)])
print(indep, rels)

# Intersection of subspaces (given by bases) in the same space
U = [f_alpha, f_beta]
V = [U2g, M11.degeneracy_map(22,1)(T2f)]
W = Intersection(U, V)
print("Intersection basis size:", len(W))

# Containment
contained, coeffs = ContainedIn(f_alpha, [f_beta, U2g])
print("Contained?", contained)
if contained: print("Coeffs:", coeffs)
"""
