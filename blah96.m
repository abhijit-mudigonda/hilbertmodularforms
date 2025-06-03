load "config.m";

SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

N := 3*ZF;
k := [2,2,2];
M := GradedRingOfHMFs(F, 40);
Mk := HMFSpace(M, N, k);
CuspFormBasis(Mk);
