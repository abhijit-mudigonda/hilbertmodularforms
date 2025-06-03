load "config.m";

SetVerbose
// Test whether the square of a [2,2,3] form is in the [4,4,6] space

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

N := ideal<ZF | [43, 0, 0], [15, 1, 0]>;
k := [2,2,3];
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.1;

assert Order(chi) eq 2;
assert HeckeCharLabel(chi) eq "1.-2.-1.1_43.2_2u1u1.2.3u";

M := GradedRingOfHMFs(F, 111);
M223 := HMFSpace(M, N, k, chi);
M446 := HMFSpace(M, N, [4,4,6]);

S223 := CuspFormBasis(M223);
S446 := CuspFormbasis(M446);

S223_squared := [f^2 : f in S223];
#Intersection(S223_squared, S446);

