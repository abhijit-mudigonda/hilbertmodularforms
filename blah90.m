load "config.m";

SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

k := [2,2,2];
N := 3*ZF;

M := GradedRingOfHMFs(F, 100);
Mk := HMFSpace(M, N, k);
Sk := CuspFormBasis(Mk);


