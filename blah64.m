load "config.m";

SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);
N := 13*ZF;

k := [2,2,3];
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.3;

B := QuaternionAlgebra<F | -3, -12*F.1>;
O := MaximalOrder(B);
Gamma := FuchsianGroup(O);

// M := HeckeMatrix2(Gamma, N, "Infinity", k, chi);
// print Nrows(M);


