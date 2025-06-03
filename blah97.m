import "ModFrmHil/diamond.m" : AlgQuatEltLabel, SaveHeckeMatrix; 
import "ModFrmHil/level.m" : InducedH1;
load "config.m";
SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

B := QuaternionAlgebra<F | -3, -12*F.1>;
O := MaximalOrder(B);
Gamma := FuchsianGroup(O);

N := 12*ZF;
k := [2,2,3];
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.4;
assert Order(chi) eq 2;
assert HeckeCharLabel(chi) eq "1.-2.-1.1_1728.1_2u1.1.1.1u1.2.3u";

// pp := Factorization(83*ZF)[1][1];
// SaveHeckeMatrix(Gamma, N, pp, k, chi);

X := cIdealDatum(Gamma, N : chi:=chi);
InducedH1(X, X, [2,2,3]);
