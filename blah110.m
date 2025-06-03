load "config.m";

import "ModFrmHil/diamond.m" : update_coset_reps;

SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

B := QuaternionAlgebra<F | -3, -12*F.1>;
O := MaximalOrder(B);
Gamma := FuchsianGroup(O);
_ := Group(Gamma);
Gamma`ShimFDSidepairsDomain;
