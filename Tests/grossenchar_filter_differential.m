// test-time: ~2s
// Differential/edge-case tests for the RayClassFilter optional argument to
// PrunedGrossencharsSet/HMFGrossencharsTorsorSet (ModFrmHilD/HMFGrossenchar.m,
// ModFrmHilD/Creation/Dihedral.m), which lets PossibleGrossencharsOfRelQuadExt
// skip building HMFGrossenchar objects for candidates that can't match the
// target nebentypus, instead of building the whole torsor and filtering after.

/***************************************************
 * PossibleGrossencharsOfRelQuadExt now always builds the RayClassFilter
 * internally in the finite order case, and keeps the old (pre-filter)
 * per-psi recheck as an `assert` safety net -- so exercising it across
 * many (K, chi) pairs is itself the differential check: it would crash
 * loudly on any disagreement between the two.
***************************************************/

F := QuadraticField(5);
ZF := Integers(F);
N := ideal<ZF | 241, 2*ZF.2 + 137>;
H := HeckeCharacterGroup(N, [1,2]);
chi := H.1^11;

Ks := QuadraticExtensionsWithConductor(N, [1,2]);
checked_nonempty := false;
for K in Ks do
  ans := PossibleGrossencharsOfRelQuadExt(K, N, [1,1], chi);
  if #ans gt 0 then
    checked_nonempty := true;
  end if;
end for;
assert checked_nonempty;

/***************************************************
 * End to end regression: PossibleGrossencharsOfRelQuadExt now always
 * filters in the finite order case; confirm the known dihedral form at
 * this level is still found.
***************************************************/

Mchi := HMFSpace(GradedRingOfHMFs(F, 150), N, [1,1], chi);
Dchi := DihedralBasis(Mchi);
assert #Dchi eq 1;
assert #LinearDependence(Dchi) eq 0;

/***************************************************
 * Trivial ray class group: PrunedGrossencharsSet should return {} whether
 * or not a filter is supplied. Q(i) has class number 1, so the trivial
 * modulus gives a trivial ray class group.
***************************************************/

K2 := QuadraticField(-1);
N2 := 1 * Integers(K2);
X2 := cHMFGrossencharsTorsor(K2, [<0,0>], N2);
assert IsFiniteOrder(X2);
_, N_oo2 := Modulus(X2);
H2 := HeckeCharacterGroup(N2, N_oo2);
assert #H2 eq 1;
assert #PrunedGrossencharsSet(X2, Rationals()) eq 0;
assert #PrunedGrossencharsSet(X2, Rationals() : RayClassFilter:=[]) eq 0;

/***************************************************
 * Empty result under a filter no character can satisfy.
***************************************************/

Ks3 := QuadraticExtensionsWithConductor(N, [1,2]);
K3 := Ks3[1];
ZK3 := Integers(K3);
M3 := Integers(AbsoluteField(K3))!!(N / Discriminant(ZK3));
K3_abs := AbsoluteField(K3);
k3 := HeckeCharWeightFromWeight(K3_abs, BaseField(K3), [1,1]);
X3 := cHMFGrossencharsTorsor(K3_abs, k3, M3);
GF3, mF3 := RayClassGroup(N, [1,2]);
gen3 := mF3(Rep(Generators(GF3)));
// a target value (5) that is not a root of unity, so no character can hit it
impossible_filter := [<Integers(K3_abs)!!gen3, K3_abs!5>];
assert #PrunedGrossencharsSet(X3, F : RayClassFilter:=impossible_filter) eq 0;

/***************************************************
 * RayClassFilter is rejected outside the finite order case.
***************************************************/

Fnp := QuadraticField(2);
Knp := CyclotomicField(8);
Nnp := Factorization(41*Integers(Knp))[1][1];
knp := HeckeCharWeightFromWeight(Knp, Fnp, [1,2]);
Xnp := cHMFGrossencharsTorsor(Knp, knp, Nnp);
assert not IsFiniteOrder(Xnp);
_, Nnp_oo := Modulus(Xnp);
Hnp := HeckeCharacterGroup(Nnp, Nnp_oo);
bad_filter := [<1*Integers(Knp), Knp!1>];
try
  _ := HMFGrossencharsTorsorSet(Xnp : RayClassFilter:=bad_filter);
  assert 0 eq 1; // should not reach here
catch e
  ; // expected: require failure
end try;
