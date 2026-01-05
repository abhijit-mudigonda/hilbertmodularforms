/*
Test that the ideal coefficient path and Fourier coefficient path
for computing inclusion/degeneracy maps give the same results (up to constants).

This tests the logic in line 826 of ModFrmHilD/Components.m
*/

procedure test_old_cusp_forms_modes(MkN1, MkN2 : verbose := false)
  /*
    Test that OldCuspForms computed with ideal coefficients (default for paritious)
    and Fourier coefficients give the same result (up to constants).
  */
  if not IsParitious(Weight(MkN1)) then
    // Only test for paritious weights where there's a choice
    return;
  end if;
  
  // Compute OldCuspForms both ways
  old_ideal := OldCuspForms(MkN1, MkN2 : UseFourierCoeffs:=false);
  old_fourier := OldCuspForms(MkN1, MkN2 : UseFourierCoeffs:=true);
  
  if #old_ideal ne #old_fourier then
    if verbose then
      printf "ERROR: Different number of old forms!\n";
      printf "  Ideal mode: %o forms\n", #old_ideal;
      printf "  Fourier mode: %o forms\n", #old_fourier;
    end if;
    error "Different number of old forms!";
  end if;
  
  if #old_ideal eq 0 then
    // No old forms to compare
    if verbose then
      printf "  No old forms (dimension 0)\n";
    end if;
    return;
  end if;
  
  // Check that the spaces spanned are the same
  // The two bases should span the same space, which means:
  // #LinearDependence(old_ideal cat old_fourier) should equal #old_ideal
  combined := old_ideal cat old_fourier;
  num_deps := #LinearDependence(combined);
  
  if num_deps ne #old_ideal then
    if verbose then
      printf "ERROR: Spaces do not match!\n";
      printf "  Expected %o dependencies, got %o\n", #old_ideal, num_deps;
      printf "  Ideal mode: %o forms\n", #old_ideal;
      printf "  Fourier mode: %o forms\n", #old_fourier;
    end if;
    error "OldCuspForms modes give different spaces!";
  end if;
  
  if verbose then
    printf "  OldCuspForms modes match (%o forms)\n", #old_ideal;
  end if;
end procedure;

procedure test_inclusion_for_space(Mk, d : verbose := false)
  /*
    Given a space Mk at level N and an integral ideal d,
    test that the inclusion map from Mk to the space at level N*d
    gives the same result (up to constants) using both the ideal
    coefficient path and the Fourier coefficient path.
  */
  M := Parent(Mk);
  N := Level(Mk);
  k := Weight(Mk);
  chi := Character(Mk);
  
  // Create the target space at level N*d
  Nd := N * d;
  n := Degree(BaseField(M));
  chi_extended := Extend(chi, Nd, [1..n]);
  Mk_Nd := HMFSpace(M, Nd, k, chi_extended);
  
  // Test OldCuspForms modes
  test_old_cusp_forms_modes(Mk, Mk_Nd : verbose := verbose);
  
  // Get cusp forms at level N
  Sk := CuspFormBasis(Mk);
   
  // For each form, compute its inclusion into Mk_Nd both ways
  for f in Sk do
    // Compute inclusion using standard path (ideal coefficients for paritious)
    g_ideal := Inclusion(f, Mk_Nd, d);
    
    // Compute inclusion using forced coefficient path
    g_coeff := Inclusion(f, Mk_Nd, d : ForceCoeffPath:=true);
    
    // Check that they are the same up to a constant
    // i.e., they should be linearly dependent
    num_deps := #LinearDependence([g_ideal, g_coeff]);
    
    if num_deps ne 1 then
      if verbose then
        printf "ERROR: Inclusion paths do not match!\n";
        printf "  Field: %o\n", BaseField(M);
        printf "  Level: %o\n", N;
        printf "  Weight: %o\n", k;
        printf "  Nebentypus: %o\n", chi;
        printf "  Ideal d: %o\n", d;
        printf "  Number of linear dependencies: %o (expected 1)\n", num_deps;
      end if;
      error "Inclusion paths do not match!";
    end if;
  end for;
end procedure;

procedure test_inclusion_paths(F, N, k, test_ideals : chi:=false, verbose:=false)
  prec := 150;
  M := GradedRingOfHMFs(F, prec);
  
  if chi cmpeq false then
    Mk := HMFSpace(M, N, k);
  else
    Mk := HMFSpace(M, N, k, chi);
  end if;
  
  for d in test_ideals do
    test_inclusion_for_space(Mk, d : verbose := verbose);
  end for;

  print "-----------------------------------------------------------------------";
end procedure;

// Test cases

F := QuadraticField(5);
ZF := Integers(F);

k := [2,2];
N := Factorization(7*ZF)[1][1];
test_ideals := [2*ZF, 7*ZF];
test_inclusion_paths(F, N, k, test_ideals : verbose := true);

N := ideal<ZF | 11, 2*ZF.2 + 3>;
k := [2,4];
test_ideals := [2*ZF, N];
test_inclusion_paths(F, N, k, test_ideals : verbose := true);

/*
N := 6*ZF;
k := [3,3];
H := HeckeCharacterGroup(N, [1,2]);
chi := H.1;
assert Order(chi) eq 2;
test_ideals := [Factorization(5*ZF)[1][1], 2*ZF];
test_inclusion_paths(F, N, k, test_ideals : chi:=chi, verbose:=true);
*/

printf "All tests passed!\n";


