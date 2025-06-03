import "ModFrmHil/diamond.m" : SaveHeckeMatrix, GetHeckeMatrix;
load "config.m";

SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

N := 5*ZF;
k := [2,2,4];
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.0;

/*
N := ideal<ZF | [43, 0, 0], [15, 1, 0]>;
k := [2,2,3];
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.1;
assert Order(chi) eq 2;
assert HeckeCharLabel(chi) eq "1.-2.-1.1_43.2_2u1u1.2.3u";
*/

M := HilbertCuspForms(F, N, chi, k);

pp := 2*ZF;

// hm := GetHeckeMatrix(M, pp : SaveAndLoad:=true);
// SaveHeckeMatrix(M, pp);
hm_2, pi_2 := GetHeckeMatrix(M, pp : SaveAndLoad:=true);
// hm_3 := GetHeckeMatrix(M, pp : SaveAndLoad:=true);
hm_2;
pi_2;

