load "config.m";

import "ModFrmHil/diamond.m" : GetHeckeMatrix;
SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

B := QuaternionAlgebra<F | -3, -12*F.1>;
O := MaximalOrder(B);
Gamma := FuchsianGroup(O);

N := ideal<ZF | [56, 0, 0], [16, 8, 0]>;
k := [2,2,3];
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.5;
assert Order(chi) eq 2;
assert HeckeCharLabel(chi) eq "1.-2.-1.1_3584.1_2u1.1.0.0.0u1.2.3u";

/*

M := HilbertCuspForms(F, N, chi, k);

Tinf := GetHeckeMatrix(M, "Infinity" : shapiro_trick:=true);
*/



