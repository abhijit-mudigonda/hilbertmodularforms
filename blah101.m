load "config.m";

import "ModFrmHil/level.m" : CompleteRelationFromUnit;
SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

B := QuaternionAlgebra<F | -3, -12*F.1>;
O := MaximalOrder(B);
Gamma := FuchsianGroup(O);
U, mm, m := Group(Gamma);
X := cIdealDatum(Gamma, 6*ZF);

O3 := Order(O, 3*ZF);
Gamma3 := FuchsianGroup(O3);
pp := Factorization(7*ZF)[1][1];
U3, mm3, m3 := Group(Gamma);
X3 := cIdealDatum(Gamma3, pp);


