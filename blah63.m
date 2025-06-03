load "config.m";

SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

epses := UnitsGenerators(F);
eps := epses[1];
assert [Sign(Evaluate(eps, v)) : v in InfinitePlaces(F)] eq [1, 1, -1];
S<y> := PolynomialRing(F);

K := ext<F | y^2 - eps>;
K2 := ext<F | y^2 - eps>;


V := VectorSpace(K, 2);
V2 := VectorSpace(K2, 2);

W := sub<V | V.1 + 3*V.2>;

W2 := sub<V2 | [V2!w : w in Basis(W)]>;


