load "config.m";

SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

B := QuaternionAlgebra<F | -3, -12*F.1>;
O := MaximalOrder(B);
print Basis(O);

O3 := Order(O, 3*ZF);
O3_basis := Basis(O3);
print O3_basis;
O3_new := Order(O, O3_basis);


