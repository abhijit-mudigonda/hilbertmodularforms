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

pp := Factorization(7*ZF)[1][1];
// pp := 2*ZF;

M := HeckeMatrix2(Gamma, N, pp, k, chi);


O_3 := Order(O, 3*ZF);
Gamma_3 := FuchsianGroup(O_3);
N_3 := 4*ZF;
H_3 := HeckeCharacterGroup(N_3, [1,2,3]);
chi_3 := Restrict(chi, H_3);

M_3 := HeckeMatrix2(Gamma_3, N_3, pp, k, chi_3);
