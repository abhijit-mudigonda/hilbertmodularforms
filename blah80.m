load "config.m";

SetVerbose("ModFrmHil", 3);
SetVerbose("HilbertModularForms", 3);

F := QuadraticField(5);
ZF := Integers(F);
N := 7*ZF;
k := [2,4];

M := GradedRingOfHMFs(F, 144);
Mk := HMFSpace(M, N, k);
Sk := NewCuspFormBasis(Mk);
