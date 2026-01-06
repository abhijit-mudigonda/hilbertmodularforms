from sage.all import *

"""
Test if Atkin-Lehner operator W_p interchanges degeneracy maps alpha and beta.

For each Galois orbit U of newforms at level N, checks if:
    W_p(alpha(U, p)) = beta(U, p)
    W_p(beta(U, p)) = alpha(U, p)

where alpha and beta are degeneracy maps from level N to level Np.
"""

def dirichlet_character_from_lmfdb_label(label):
    """
    Get a Dirichlet character from its LMFDB label.
    
    LMFDB labels have the format: N.n where N is the modulus and n is the index.
    The index n corresponds to the Conrey label, which labels characters by
    the smallest positive integer in their Galois orbit.
    
    Args:
        label: string in format "N.n" (e.g., "26.2")
    
    Returns:
        Dirichlet character
    
    Examples:
        >>> chi = dirichlet_character_from_lmfdb_label("26.2")
        >>> chi.modulus()
        26
        >>> chi.order()
        2
    """
    parts = label.split('.')
    if len(parts) != 2:
        raise ValueError(f"Invalid LMFDB label format: {label}. Expected 'N.n'")
    
    try:
        N = Integer(parts[0])
        n = Integer(parts[1])
    except:
        raise ValueError(f"Invalid LMFDB label: {label}. Could not parse as integers.")
    
    # In Sage, DirichletGroup uses a different indexing than LMFDB's Conrey labels
    # LMFDB uses Conrey labels: the character that sends the generator to zeta^n
    # We need to find the character in Sage's DirichletGroup that corresponds to this
    
    # For the Conrey labeling, character number n sends the standard generator to zeta^n
    # where zeta is a primitive root mod N
    
    # LMFDB uses Conrey labels: for modulus N, the Conrey character m is defined by
    # chi_m(a) = exp(2*pi*i * ind_g(a) * ind_g(m) / phi(N))
    # where g is the smallest primitive root mod N and ind_g is the discrete log base g
    
    # In Sage, we can find the character by its Conrey number
    G = DirichletGroup(N)
    
    # Try using from_conrey if available (Sage 9.5+)
    if hasattr(G, 'from_conrey'):
        return G.from_conrey(n)
    
    # Fallback: find the character by checking conrey_number() for each character
    for chi in G:
        if hasattr(chi, 'conrey_number') and chi.conrey_number() == n:
            return chi
    
    raise ValueError(f"Could not find character with Conrey label {n} at modulus {N}")

def test_atkin_lehner_degeneracy(N, p, k=2, character=None, sign=0, verbose=True):
    """
    Test if W_p interchanges alpha(U) and beta(U) for all Galois orbits U at level N.
    
    Args:
        N: level of the starting space (p may or may not divide N)
        p: prime for degeneracy maps and Atkin-Lehner operator
        k: weight (default 2)
        character: Dirichlet character (default None = trivial character)
        sign: sign for modular symbols (0, 1, or -1)
        verbose: print detailed information (if False, just returns boolean)
    
    Returns:
        If verbose=True: Dictionary with 'all_pass', 'orbits', 'summary'
        If verbose=False: Boolean (True if all orbits pass, False otherwise)
    """
    if verbose:
        print(f"{'=' * 80}")
        char_str = "trivial" if character is None else str(character)
        print(f"Testing: Level {N} -> {N*p}, Weight {k}, Character {char_str}, Prime p = {p}")
        print(f"{'=' * 80}")
    
    # Check if p divides N
    p_divides_N = (N % p == 0)
    if verbose:
        if p_divides_N:
            print(f"Note: p = {p} divides N = {N}")
        else:
            print(f"Note: p = {p} does not divide N = {N}")
        print()
    
    # Determine the character to use
    if character is None:
        # Trivial character
        MS_N = ModularSymbols(N, k, sign=sign).cuspidal_submodule()
        MS_Np = ModularSymbols(N*p, k, sign=sign).cuspidal_submodule()
    else:
        # Non-trivial character
        # For degeneracy maps, we need to extend the character from level N to level Np
        MS_N = ModularSymbols(character, k, sign=sign).cuspidal_submodule()
        # Extend character to level Np by composing with the natural map Z/NpZ -> Z/NZ
        MS_Np = ModularSymbols(character.extend(N*p), k, sign=sign).cuspidal_submodule()
    
    if verbose:
        print(f"Cuspidal Modular Symbols:")
        print(f"  Level {N}: dimension {MS_N.dimension()}")
        print(f"  Level {N*p}: dimension {MS_Np.dimension()}")
    
    # Get newforms at level N
    newforms_N = MS_N.new_subspace()
    if verbose:
        print(f"\nNew subspace at level {N}: dimension {newforms_N.dimension()}")
    
    # Decompose into Galois orbits
    decomp = newforms_N.decomposition()
    
    if verbose:
        print(f"Number of Galois orbits: {len(decomp)}")
        if len(decomp) == 0:
            print("(No newforms at this level)")
        print()
    
    if len(decomp) == 0:
        # No newforms, so vacuously true
        if verbose:
            print(f"{'=' * 80}")
            print("RESULT: ✓ (vacuously true - no newforms)")
            print(f"{'=' * 80}\n")
            return {
                'all_pass': True,
                'orbits': [],
                'summary': 'No newforms at this level'
            }
        else:
            return True
    
    # Degeneracy maps: alpha (d=1) and beta (d=p)
    alpha_map = MS_N.degeneracy_map(N*p, 1)
    beta_map = MS_N.degeneracy_map(N*p, p)
    
    # Atkin-Lehner operator W_p at level Np
    try:
        W_p = MS_Np.atkin_lehner_operator(p)
    except (NotImplementedError, ValueError) as e:
        if verbose:
            print(f"ERROR: Atkin-Lehner operator W_{p} not available at level {N*p}")
            print(f"Reason: {e}")
            print(f"{'=' * 80}\n")
            return {
                'all_pass': None,
                'orbits': [],
                'summary': f'W_{p} not available: {e}'
            }
        else:
            return None  # Return None to indicate "not applicable" rather than False
    
    orbit_results = []
    
    for i, U in enumerate(decomp):
        if verbose:
            print(f"{'-' * 80}")
            print(f"Galois Orbit {i+1}:")
            print(f"  Dimension: {U.dimension()}")
        
        # Get basis of U
        basis_U = U.basis()
        
        # Apply alpha and beta to each basis element
        alpha_U_gens = [alpha_map(f) for f in basis_U]
        beta_U_gens = [beta_map(f) for f in basis_U]
        
        # Apply W_p to the images
        W_alpha_U_gens = [W_p(f) for f in alpha_U_gens]
        W_beta_U_gens = [W_p(f) for f in beta_U_gens]
        
        # Check if span(W_p(alpha(U))) = span(beta(U)) and span(W_p(beta(U))) = span(alpha(U))
        try:
            # Express each generator as coordinates in the ambient space
            def get_coords(element):
                return element.element()
            
            # Build coordinate matrices
            W_alpha_coords = matrix([get_coords(f) for f in W_alpha_U_gens])
            beta_coords = matrix([get_coords(f) for f in beta_U_gens])
            W_beta_coords = matrix([get_coords(f) for f in W_beta_U_gens])
            alpha_coords = matrix([get_coords(f) for f in alpha_U_gens])
            
            # Check if row spaces are equal
            relation1 = (W_alpha_coords.row_space() == beta_coords.row_space())
            relation2 = (W_beta_coords.row_space() == alpha_coords.row_space())
            
        except Exception as e:
            if verbose:
                print(f"  Error checking relations: {e}")
            relation1 = False
            relation2 = False
        
        both_hold = relation1 and relation2
        
        if verbose:
            print(f"\n  Checking relations:")
            print(f"    W_{p}(alpha(U)) == beta(U):  {relation1}")
            print(f"    W_{p}(beta(U)) == alpha(U):  {relation2}")
            if both_hold:
                print(f"    ✓ Both relations hold!")
            else:
                print(f"    ✗ Relations do NOT hold")
            print()
        
        orbit_results.append({
            'dimension': U.dimension(),
            'W_p_alpha_equals_beta': relation1,
            'W_p_beta_equals_alpha': relation2,
            'both_hold': both_hold
        })
    
    all_pass = all(r['both_hold'] for r in orbit_results)
    
    if verbose:
        print(f"{'=' * 80}")
        print("RESULT:")
        if all_pass:
            print(f"✓ All {len(orbit_results)} Galois orbit(s): W_{p} interchanges alpha and beta")
        else:
            num_pass = sum(1 for r in orbit_results if r['both_hold'])
            print(f"✗ Only {num_pass}/{len(orbit_results)} Galois orbit(s) satisfy the relation")
        print(f"{'=' * 80}\n")
    
    summary = f"Level {N} -> {N*p} (wt {k}): {'PASS' if all_pass else 'FAIL'}"
    
    if verbose:
        return {
            'all_pass': all_pass,
            'orbits': orbit_results,
            'summary': summary
        }
    else:
        return all_pass


# Example 1: Trivial character
print("Example 1: Trivial character at level 11, weight 4")
test_atkin_lehner_degeneracy(N=11, p=2, k=4, sign=0, verbose=True)

print("\n\n")

# Example 1: Trivial character

print("Example 1.1: Trivial character at level 22, weight 2")
test_atkin_lehner_degeneracy(N=11, p=2, k=2, sign=0, verbose=True)
print("\n\n")

print("Example 1.2: Trivial character at level 22, weight 2")
test_atkin_lehner_degeneracy(N=75, p=5, k=2, sign=0, verbose=True)
print("\n\n")

# Example 2: Non-trivial character at level 26
# The quadratic character at level 26 has Conrey label 25
print("Example 2: Quadratic character at level 26, weight 2 (Conrey label 26.25)")
chi = dirichlet_character_from_lmfdb_label("26.25")
print(f"Using character: {chi}, order {chi.order()}")
print()
test_atkin_lehner_degeneracy(N=26, p=2, k=2, character=chi, sign=0, verbose=True)

print("\n\n")

# Example 3: Quartic character at level 16
print("Example 3: Quartic character at level 16, weight 2 (Conrey label 16.5)")
chi = dirichlet_character_from_lmfdb_label("16.5")
print(f"Using character: {chi}, order {chi.order()}")
print()
test_atkin_lehner_degeneracy(N=16, p=2, k=2, character=chi, sign=0, verbose=True)

print("\n\n")

# Example 4: Test verbose=False mode
print("Example 4: Testing verbose=False mode")
result = test_atkin_lehner_degeneracy(N=7, p=2, k=4, sign=0, verbose=False)
print(f"Result (boolean): {result}")
