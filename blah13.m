load "config.m";

F := QuadraticField(2);
ZF := Integers(F);
M := GradedRingOfHMFs(F, 200);
k := [2, 3];

N := ideal<ZF | 28, 2*ZF.2 - 8>;
H := HeckeCharacterGroup(N, [1,2]);
chi := H.2;
assert Order(chi) eq 2;
assert IsCompatibleWeight(chi, k);

Mk := HMFSpace(M, N, k, chi);
Sk := CuspFormBasis(Mk : StableOnly:=true);

