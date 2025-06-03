load "config.m";

import "ModFrmHil/level.m" : CompleteRelationFromUnit;
SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

B := QuaternionAlgebra<F | -3, -12*F.1>;
O := MaximalOrder(B);
Gamma := FuchsianGroup(O);


pp := Factorization(7*ZF)[1][1];
// pp := Factorization(13*ZF)[1][1];
// pp := 3*ZF;
// pp := 2*ZF;
k := [2,2,4];
// k := [2,2,2];

O2 := Order(O, 3*ZF);
Gamma2 := FuchsianGroup(O2);
N2 := 2*ZF;
// N2 := 1*ZF;
chi2 := HeckeCharacterGroup(N2, [1,2,3]).0;
// T2 := HeckeMatrix2(Gamma2, N2, pp, k, chi2);

X := cIdealDatum(Gamma, 6*ZF);
X2 := cIdealDatum(Gamma2, 2*ZF);

U, _, m := Group(Gamma);
U2, _, m2 := Group(Gamma2);
CompleteRelationFromUnit(X, m(U.2), k);
/*
N := 6*ZF;
chi := HeckeCharacterGroup(N, [1,2,3]).0;
T := HeckeMatrix2(Gamma, N, pp, k, chi);

CharacteristicPolynomial(T2) eq CharacteristicPolynomial(T);
*/

