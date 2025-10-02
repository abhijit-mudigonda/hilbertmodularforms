// tests [1,2] computation at a single level, character
// for debug

load "config.m";

F := QuadraticField(2);
ZF := Integers(F);
M := GradedRingOfHMFs(F, 500);
k := [1, 2];

N := ideal<ZF | 28, 2*ZF.2 - 8>;
H := HeckeCharacterGroup(N, [1,2]);
chi := H.1;
assert Order(chi) eq 2;
assert HeckeCharLabel(chi) eq "-2.0.1_56.1_2u0.1u1.2u";
assert IsCompatibleWeight(chi, k);

Mk := HMFSpace(M, N, k, chi);
Sk := CuspFormBasis(Mk : StableOnly:=true);
