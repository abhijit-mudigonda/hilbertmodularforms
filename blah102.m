load "config.m";

import "ModFrmHil/diamond.m" : GetHeckeMatrix;
SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

B := QuaternionAlgebra<F | -3, -12*F.1>;
O := MaximalOrder(B);
Gamma := FuchsianGroup(O);

/*
N := 3*ZF;
k := [2,2,4];
chi := HeckeCharacterGroup(N, [1,2,3]).0;
*/

N := ideal<ZF | [301, 0, 0], [58, 1, 0] >;
k := [2,2,3];
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.2;
assert Order(chi) eq 2;
assert HeckeCharLabel(chi) eq "1.-2.-1.1_301.2_2u1.1u1.2.3u";

M := HilbertCuspForms(F, N, chi, k);

Tinf := GetHeckeMatrix(M, "Infinity" : shapiro_trick:=true);
T := AssociativeArray();
for pp in [Factorization(7*ZF)[1][1], 2*ZF, Factorization(13*ZF)[1][1], 3*ZF] do
  T[Norm(pp)] := GetHeckeMatrix(M, pp : shapiro_trick:=true);
end for;

assert IsZero(T[7]*T[13]-T[13]*T[7]);
assert IsZero(T[8]*T[13]-T[13]*T[8]);
assert IsZero(T[8]*T[27]-T[27]*T[8]);
assert IsZero(T[27]*T[13]-T[13]*T[27]);

print T[7];
print T[8];

/*
O_3 := Order(O, 3*ZF);
Gamma_3 := FuchsianGroup(O_3);
H_3 := HeckeCharacterGroup(2*ZF, [1,2,3]);
chi_3 := Restrict(chi, H_3);
U_m, mm_m, m_m := Group(Gamma_3);

pp := Factorization(7*ZF)[1][1];
T7 := HeckeMatrix2(Gamma_3, 2*ZF, pp, k, chi_3);

pp := Factorization(13*ZF)[1][1];
T13 := HeckeMatrix2(Gamma_3, 2*ZF, pp, k, chi_3);

pp := 2*ZF;
T2 := HeckeMatrix2(Gamma_3, 2*ZF, pp, k, chi_3);

Tinf := HeckeMatrix2(Gamma_3, 2*ZF, "Infinity", k, chi_3);

pp := 3*ZF;
T3 := HeckeMatrix2(Gamma, 6*ZF, pp, k, chi);
*/
