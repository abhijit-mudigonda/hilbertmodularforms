load "config.m";

import !"Geometry/ModFrmHil/indefinite.m" : ElementOfNormMinusOne;

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

B := QuaternionAlgebra<F | -3, -12*F.1>;
O := MaximalOrder(B);
Gamma := FuchsianGroup(O);

mu := ElementOfNormMinusOne(O);
eps := Norm(mu);
auts := AutsOfKReppingEmbeddingsOfF(F, F);
k := [1,1,2];
k0 := Max(k);
// TODO abhijitm
// There are several assumptions being made here
// - the quaternion algebra is split at exactly one infinite prime
// - at this prime, the weight is k0

// This value is the the value of the determinant character of the 
// weight representation on the chosen "ElementOfNormMinusOne" 
// (which represents complex conjugation). The point is that the + and -
// subspaces are not defined 
z := &*[auts[i](eps)^(k0 - k[i]) : i in [1 .. #auts]];

// _, K, _ := Splittings(B);
R<x> := PolynomialRing(F);
poly := x^2 - z;
if IsIrreducible(poly) then
  // TODO abhijitm I'm a little bit worried about returning a 
  // relative extension here instead of an absolute one...
  K := ext<F | x^2 - z>;
end if;

