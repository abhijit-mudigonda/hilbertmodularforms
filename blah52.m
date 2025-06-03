import !"Geometry/ModFrmHil/indefinite.m" : ElementOfNormMinusOne;

load "config.m";

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

N := 4*ZF;
// N := 8*ZF;
k := [2,2,3];
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.3;
// chi := H.4;
assert Order(chi) eq 2;
assert HeckeCharLabel(chi) eq "1.-2.-1.1_64.1_2u1.1.1u1.2.3u";
// assert HeckeCharLabel(chi) eq "1.-2.-1.1_512.1_2u1.1.1.0u1.2.3u";

B := QuaternionAlgebra<F | -3, -12*F.1>;
// B := QuaternionAlgebra<F | -F.1^2 + F.1 - 1, -8*F.1^2 + 4*F.1 + 16>;
O := MaximalOrder(B);
Gamma := FuchsianGroup(O);

hm_inf := HeckeMatrix2(Gamma, N, "Infinity", k, chi);

/*
K := BaseRing(hm_inf);
chi_inf := CharacteristicPolynomial(hm_inf);
L := ext<K | Factorization(chi_inf)[1][1]>;

mu := ElementOfNormMinusOne(O);
eps := Norm(mu);
vs := InfinitePlaces(F);
[Evaluate(eps, v) : v in vs];

eps_1 := Evaluate(eps, vs[1]);
eps_2 := Evaluate(eps, vs[2]);

Sqrt(eps_1) * Sqrt(eps_2);

eps_1 * eps_2;

[Real(Evaluate(L.1, w)) : w in InfinitePlaces(L)];

[Real(Evaluate(L.1^2, w)) : w in InfinitePlaces(L)];
*/


/*
L := NumberField(x^48 - 12*x^46 + 30*x^44 + 212*x^42 - 1005*x^40 - 108*x^38 +10056*x^36 - 16476*x^34 - 16659*x^32 + 101124*x^30 - 95058*x^28 - 146832*x^26 +365175*x^24 - 146832*x^22 - 95058*x^20 + 101124*x^18 - 16659*x^16 - 16476*x^14 +10056*x^12 - 108*x^10 - 1005*x^8 + 212*x^6 + 30*x^4 - 12*x^2 + 1);
K := BaseRing(hm);
IsSubfield(K, L);
hm_L := ChangeRing(hm, L);

[Evaluate(Determinant(hm), v) : v in InfinitePlaces(K)];
*/

hecke_mtrxs := AssociativeArray();

for pp in PrimesUpTo(20, F) do
  print "pp", IdealOneLine(pp);
  hecke_mtrx := HeckeMatrix2(Gamma, N, pp, k, chi);
  hecke_mtrxs[pp] := hecke_mtrx;
  print Eigenvalues(hecke_mtrx);
end for;

/*
function g(M, aut)
  A := [[aut(M[i][j]) : j in [1 .. Ncols(M)]] : i in [1 .. Nrows(M)]];
  return Matrix(A);
end function;

auts := AutsOfKReppingEmbeddingsOfF(F, K);

blorp := &*[g(hm_inf, aut) : aut in auts];
*/
