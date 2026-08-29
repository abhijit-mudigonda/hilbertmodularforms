/***************************************************
* Regression tests for the tight per-divisor RayClassFilter search
* (TightRayClassFilteredGrossencharsSet, PruneConjugatePairs), which
* PossibleGrossencharsOfRelQuadExt's finite-order branch (e.g. weight
* [1,1]) uses instead of materializing HeckeCharacterGroup(Modulus(X))
* in full.
***************************************************/

/***************************************************
* Q(sqrt5), N=23, weight [1,1], chi = H.1^11 (Tests/wt1lvl23.m's exact
* setup): golden value cross-checked directly against the pre-rewrite
* full-enumeration implementation.
***************************************************/
F := QuadraticField(5);
M := GradedRingOfHMFs(F, 100);
ZF := Integers(F);
N := 23*ZF;
H := HeckeCharacterGroup(N, [1,2]);
chi := H.1^11;
Mk := HMFSpace(M, N, [1,1], chi);
Dk := DihedralBasis(Mk);
assert #Dk eq 1;

/***************************************************
* Q(sqrt6), N=56/57, trivial chi, weight [1,1]: fix(Dihedral): factor
* primes in the absolute order, not the relative order. K's relative
* defining polynomial is non-integral here, which is exactly the
* condition under which the old relative-order route (still present in
* MaximalIdealsOfNormDividing's counterpart before this rewrite) could
* throw "Ideal is fractional" or silently return a wrong ideal. Both
* levels have no weight-[1,1] dihedral forms with trivial character at
* any auxiliary field -- the assertion is that this completes cleanly
* and stably, not on a specific nonzero count.
***************************************************/
F6 := QuadraticField(6);
ZF6 := Integers(F6);
for Nval in [56, 57] do
  N6 := Nval*ZF6;
  chi6 := HeckeCharacterGroup(N6, [1,2])!1;
  total := 0;
  for K in QuadraticExtensionsWithConductor(N6, [1,2]) do
    total +:= #PossibleGrossencharsOfRelQuadExt(K, N6, [1,1], chi6);
  end for;
  assert total eq 0;
end for;
