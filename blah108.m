load "config.m";

import "ModFrmHil/diamond.m" : GetHeckeMatrix, update_coset_reps, shapiro_matrix;
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

qq := Factorization(7*ZF)[1][1];
O_qq := Order(O, qq);
Gamma_qq := FuchsianGroup(O_qq);
U_m, mm_m, m_m := Group(Gamma_qq);


H_qq := HeckeCharacterGroup(N / qq, [1,2,3]);
chi_qq := Restrict(chi, H_qq);

X := cIdealDatum(Gamma, N : chi:=chi);
X_qq := cIdealDatum(Gamma_qq, N / qq : chi:=chi_qq, residue_map:=X`ResidueMap);
subgp_coset_idxs := update_coset_reps(X, X_qq);

Tinf := HeckeMatrix2(Gamma, N, "Infinity", k, chi);
print "noice -----------------";
Tinf_qq := HeckeMatrix2(Gamma_qq, N / qq, "Infinity", k, chi_qq);

pp := 2*ZF;
T2 := HeckeMatrix2(Gamma, N, pp, k, chi);
T2_qq := HeckeMatrix2(Gamma_qq, N / qq, pp, k, chi_qq);

pp := Factorization(7*ZF)[1][1];
T7 := HeckeMatrix2(Gamma, N, pp, k, chi);
// T7_qq := HeckeMatrix2(Gamma_qq, N / qq, pp, k, chi_qq);
S := shapiro_matrix(X, X_qq, k, subgp_coset_idxs);
