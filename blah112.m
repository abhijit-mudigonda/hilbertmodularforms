load "config.m";

SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

M := GradedRingOfHMFs(F, 200);
N := Factorization(13*ZF)[1][1];
k := [4,4,4];
Mk := HMFSpace(M, N, k);
Sk := CuspFormBasis(Mk);
Ek := Eigenbasis(Mk, Sk : P:=10);
