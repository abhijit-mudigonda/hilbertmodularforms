load "config.m";

import "ModFrmHil/definite.m" : HeckeOperatorDefiniteBig;
import !"Geometry/ModFrmHil/hecke.m" : restriction;

MAX_PRIME := 50;

SetSeed(1337);
F := QuadraticField(2);
ZF := Integers(F);

k := [2, 3];
N := LMFDBIdeal(F, "32.1");
H := HeckeCharacterGroup(N, [1,2]);
chi := H.3;
print HeckeCharLabel(chi);
assert Order(chi) eq 2;

// COMPUTE USING DEFINITE METHOD

print "-------- DEFINITE -----------";

B_def := QuaternionAlgebra(1*ZF, InfinitePlaces(F) : Optimized);
O_def := MaximalOrder(B_def);

// def_hm_inf := HeckeOperatorDefiniteBig(M, "Infinity");

M := HilbertCuspForms(F, N, chi, k);

def_hecke_mtrxs := AssociativeArray();
for pp in PrimesUpTo(MAX_PRIME, F) do
  print "--------";
  print "pp", Norm(pp), IdealOneLine(pp);
  hecke_mtrx, p_rep := HeckeOperatorDefiniteBig(M, pp);
  print "p_rep", p_rep;
  res_hecke_mtrx := restriction(M, hecke_mtrx);
  def_hecke_mtrxs[pp] := <res_hecke_mtrx, p_rep>;
  print CharacteristicPolynomial(res_hecke_mtrx);
end for;

// COMPUTE USING HECKE STABILITY

/*
print "------------ HECKE STABILITY ------------";
GRing := GradedRingOfHMFs(F, 500);
Mk := HMFSpace(GRing, N, k, chi);
Sk := HeckeStabilityCuspBasis(Mk : prove:=false, stable_only:=true);
// in the [3,3] case I think
// assert #Sk eq 4;
// Ek := EisensteinBasis(Mk);
// assert #Ek eq 4;

for pp in primes do
  print "--------";
  print <Norm(pp), IdealOneLine(pp)>;
  print "definite";
  print Factorization(CharacteristicPolynomial(def_hecke_mtrxs[pp]));
  print "hecke stability";
  print Factorization(CharacteristicPolynomial(HeckeMatrix(Sk, pp)));
  print CharacteristicPolynomial(def_hecke_mtrxs[pp]) eq CharacteristicPolynomial(HeckeMatrix(Sk, pp));
end for;
*/
