load "config.m";

SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

prec := 100;
M := GradedRingOfHMFs(F, prec);
N := 3*ZF;
k := [2,2,4];

Mk := HMFSpace(M, N, k);
Sk := CuspFormBasis(Mk);
