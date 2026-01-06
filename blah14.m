load "config.m";
verbose := true;

F := QuadraticField(5);
ZF := Integers(F);
prec := 150;
M := GradedRingOfHMFs(F, prec);

N := ideal<ZF | 11, 2*ZF.2 + 3>;
k := [2,4];
test_ideals := [2*ZF, N];

MkN1 := HMFSpace(M, N, k);
MkN2 := HMFSpace(M, N * 2, k);

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

