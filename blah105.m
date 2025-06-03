load "config.m";

SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

M := GradedRingOfHMFs(F, 40);
N := 6*ZF;
k := [2,2,4];

Mk := HMFSpace(M, N, k);
Sk := CuspFormBasis(Mk);
Sk;

