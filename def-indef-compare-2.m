load "config.m";

import "ModFrmHil/definite.m" : HeckeOperatorDefiniteBig;
import !"Geometry/ModFrmHil/hecke.m" : restriction;

MAX_PRIME := 50;

F := QuadraticField(5);
ZF := Integers(F);

// weight [2, 2] space for sanity checking

k := [4, 4];
N := 6*ZF;
qq := N;
H := HeckeCharacterGroup(N, [1, 2]);
chi := H.0;
assert Order(chi) eq 1;


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
  print <Norm(pp), CharacteristicPolynomial(res_hecke_mtrx)>;
end for;

// COMPUTE USING INDEFINITE METHOD (away from qq = 2*ZF)

/*
print "-------- INDEFINITE -----------";
B_indef := QuaternionAlgebra(qq, [InfinitePlaces(F)[1]]);
O_indef := MaximalOrder(B_indef);
Gamma := FuchsianGroup(O_indef);

// indef_hm_inf := HeckeMatrix2(Gamma, N, "Infinity", k, chi);

indef_hecke_mtrxs := AssociativeArray();

N_indef := N / qq;

// chi_indef := HeckeCharacterGroup(N / qq, [1, 2]).;

chi_indef := HeckeCharacterGroup(N / qq, [1, 2]).1;


primes := [pp : pp in PrimesUpTo(MAX_PRIME, F) | pp ne qq];
for pp in primes do
  hecke_mtrx := HeckeMatrix2(Gamma, N_indef, pp, k, chi_indef);
  indef_hecke_mtrxs[pp] := hecke_mtrx;
  print <Norm(pp), IdealOneLine(pp), SquareRoot(CharacteristicPolynomial(hecke_mtrx))>;
end for;
*/
