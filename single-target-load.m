load "config.m";

SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

N := ideal<ZF | [43, 0, 0], [15, 1, 0]>;
k := [2,2,3];
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.1;
assert Order(chi) eq 2;
assert HeckeCharLabel(chi) eq "1.-2.-1.1_43.2_2u1u1.2.3u";

/*
N := 3*ZF;
k := [2,2,4];
chi := HeckeCharacterGroup(N, [1,2,3]).0;
*/

/*
N := 3*ZF;
k := [2,2,2];
chi := HeckeCharacterGroup(N, [1,2,3]).0;
*/

M := GradedRingOfHMFs(F, 50);
Mk := HMFSpace(M, N, k, chi);
Sk := NewCuspFormBasis(Mk);

Sk[1];

/*
M446 := HMFSpace(M, N, [4,4,6]);
S446 := CuspFormBasis(M446);

#Intersection([f^2 : f in Sk], S446);
*/
