load "config.m";

import "ModFrmHil/definite.m" : HeckeOperatorDefiniteBig;
import !"Geometry/ModFrmHil/hecke.m" : restriction;

MAX_PRIME := 50;

F := QuadraticField(5);
ZF := Integers(F);

// weight [2, 2] space for sanity checking

/*
k := [2, 2];
N := Factorization(7*ZF)[1][1];
qq := N;
H := HeckeCharacterGroup(N, [1, 2]);
chi := H.0;
assert Order(chi) eq 1;
*/

/*
// weight [3, 3] space with quadratic nebentypus

k := [3, 3];
N := 6*ZF;

// N := ideal<ZF | 4*ZF.2 - 2>;
qq := 2*ZF;
H := HeckeCharacterGroup(N, [1, 2]);
chi := H.1;
// assert HeckeCharLabel(chi) eq "-5.0.1_20.1_2u1u1.2u";
assert Order(chi) eq 2;
*/

// weight [2, 2] space with cubic nebentypus

k := [2, 2];
N := 14*ZF;
H := HeckeCharacterGroup(N, [1, 2]);
chi := H.1^2; 
assert Order(chi) eq 3;
qq := 2*ZF;

// COMPUTE USING DEFINITE METHOD

print "-------- DEFINITE -----------";

B_def := QuaternionAlgebra(1*ZF, InfinitePlaces(F) : Optimized);
O_def := MaximalOrder(B_def);

// def_hm_inf := HeckeOperatorDefiniteBig(M, "Infinity");

M := HilbertCuspForms(F, N, chi, k);

def_hecke_mtrxs := AssociativeArray();
for pp in PrimesUpTo(MAX_PRIME, F) do
  hecke_mtrx := HeckeOperatorDefiniteBig(M, pp);
  res_hecke_mtrx := restriction(M, hecke_mtrx);
  def_hecke_mtrxs[pp] := res_hecke_mtrx;
  // print <Norm(pp), IdealOneLine(pp), CharacteristicPolynomial(res_hecke_mtrx)>;
end for;

// COMPUTE USING INDEFINITE METHOD (away from qq = 2*ZF)

print "-------- INDEFINITE -----------";
B_indef := QuaternionAlgebra(qq, [InfinitePlaces(F)[1]]);
O_indef := MaximalOrder(B_indef);
Gamma := FuchsianGroup(O_indef);

_, K, _ := Splittings(B_indef);
L := Compositum(K, CyclotomicField(Order(chi)));

// indef_hm_inf := HeckeMatrix2(Gamma, N, "Infinity", k, chi);

indef_hecke_mtrxs := AssociativeArray();

N_indef := N / qq;

// chi_indef := HeckeCharacterGroup(N / qq, [1, 2]).;

chi_indef := HeckeCharacterGroup(N / qq, [1, 2]).1;


primes := [pp : pp in PrimesUpTo(MAX_PRIME, F) | pp ne qq];
for pp in primes do
  hecke_mtrx := HeckeMatrix2(Gamma, N_indef, pp, k, chi_indef : HeckeMatrixField:=L);
  indef_hecke_mtrxs[pp] := hecke_mtrx;
  // print <Norm(pp), IdealOneLine(pp), SquareRoot(CharacteristicPolynomial(hecke_mtrx))>;
end for;

// COMPUTE USING HECKE STABILITY

/*
print "------------ HECKE STABILITY ------------";
GRing := GradedRingOfHMFs(F, 500);
Mk := HMFSpace(GRing, N, k, chi);
Sk := HeckeStabilityCuspBasis(Mk : prove:=false, stable_only:=true);
// in the [3,3] case I think
assert #Sk eq 4;
// Ek := EisensteinBasis(Mk);
// assert #Ek eq 4;
*/

/* For the sanity check [2, 2] case
Sk := CuspFormBasis(Mk);
assert #Sk eq 1;
*/

for pp in primes do
  print "--------";
  print <Norm(pp), IdealOneLine(pp)>;
  print "definite";
  print Factorization(CharacteristicPolynomial(def_hecke_mtrxs[pp]));
  // assert CharacteristicPolynomial(def_hecke_mtrxs[pp]) eq CharacteristicPolynomial(HeckeMatrix(Sk, pp));
  // assert IsDivisibleBy(CharacteristicPolynomial(HeckeMatrix(Sk, pp)), SquareRoot(CharacteristicPolynomial(indef_hecke_mtrxs[pp])));
  print "indefinite";
  print Factorization(SquareRoot(CharacteristicPolynomial(indef_hecke_mtrxs[pp])));
end for;
