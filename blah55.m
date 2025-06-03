load "config.m";

SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

B := QuaternionAlgebra(1*ZF, InfinitePlaces(F)[1 .. 2]);
_, K, _ := Splittings(B);
K;
