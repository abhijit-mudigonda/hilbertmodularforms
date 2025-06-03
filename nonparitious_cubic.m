load "config.m";

SetVerbose("ModFrmHil", 3);

// Test whether the square of a [2,2,3] form is in the [4,4,6] space

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

/*
N := ideal<ZF | [301, 0, 0], [58, 1, 0] >;
k := [2,2,3];
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.2;
assert Order(chi) eq 2;
assert HeckeCharLabel(chi) eq "1.-2.-1.1_301.2_2u1.1u1.2.3u";
*/

N := 12*ZF;
k := [2,2,3];
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.4;
assert Order(chi) eq 2;
assert HeckeCharLabel(chi) eq "1.-2.-1.1_1728.1_2u1.1.1.1u1.2.3u";

assert IsCompatibleWeight(chi, k);

M := GradedRingOfHMFs(F, 100);
Mk := HMFSpace(M, N, k, chi);
Sk := CuspFormBasis(Mk);

/*
M446 := HMFSpace(M, N, [4,4,6]);
S446 := NewCuspFormBasis(M446);

assert #LinearDependence(S446) eq 0;
Sk_squared := [Sk[1]^2, Sk[1] * Sk[2], Sk[2]^2];
assert #LinearDependence(Sk_squared) eq 0;
assert #Intersection(Sk_squared, S446) eq #Sk_squared;
*/
