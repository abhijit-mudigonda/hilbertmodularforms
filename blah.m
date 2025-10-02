load "config.m";

SetVerbose("ModFrmHil", 3);
PREC := 200;
F := QuadraticField(5);
ZF := Integers(F);

N := 13*ZF;
H := HeckeCharacterGroup(N, [1,2]);
chi := H.1^6;
assert Order(chi) eq 2;

M := GradedRingOfHMFs(F, PREC);
Mk := HMFSpace(M, N, [2, 2], chi);
Sk := CuspFormBasis(Mk);

M44 := HMFSpace(M, N, [4, 4]);
S44 := CuspFormBasis(M44);
