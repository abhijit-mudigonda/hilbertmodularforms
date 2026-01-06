// 391

load "config.m";

import "ModFrmHil/definite.m" : HeckeOperatorDefiniteBig;
import !"Geometry/ModFrmHil/hecke.m" : restriction, pseudo_inverse;

MAX_PRIME := 40;

F := QuadraticField(2);
ZF := Integers(F);
k := [2,3];

chi_label := "-2.0.1_441.1_2u1.1u1.2u";
qq := ideal<ZF | 7, ZF.2 + 3>; // this won't work for indefinite stuff,
            // this is just to check definite old/new

chi := FullChiLabelToHeckeChar(chi_label);
assert Order(chi) eq 2;
assert IsCompatibleWeight(chi, k);
N := Modulus(chi);
print "Norm qq", Norm(qq);
assert N ne Conductor(chi);
assert (N / qq) subset Conductor(chi);
assert (N / (3*ZF)) subset Conductor(chi);
H := HeckeCharacterGroup(N, [1,2]);

print "-------- DEFINITE -----------";

B_def := QuaternionAlgebra(1*ZF, InfinitePlaces(F) : Optimized);
O_def := MaximalOrder(B_def);

M := HilbertCuspForms(F, N, chi, k);
M_new := NewSubspace(M, N);
_ := Dimension(M_new);  // This triggers the basis matrix computation
new_basis_matrix := M_new`basis_matrix;
new_basis_matrix_inv := M_new`basis_matrix_inv;

M_old_7 := HilbertCuspForms(F, N / qq, Restrict(chi, N / qq, [1,2]), k);
M_old_3 := HilbertCuspForms(F, N / (3*ZF), Restrict(chi, N / (3*ZF), [1,2]), k);

def_hecke_mtrxs := AssociativeArray();
def_new_hecke_mtrxs := AssociativeArray();
def_old_7_hecke_mtrxs := AssociativeArray();
def_old_3_hecke_mtrxs := AssociativeArray();

for pp in PrimesUpTo(MAX_PRIME, F) do
  hecke_mtrx := HeckeOperatorDefiniteBig(M, pp);
  def_hecke_mtrxs[pp] := hecke_mtrx;
  def_new_hecke_mtrxs[pp] := new_basis_matrix * hecke_mtrx * new_basis_matrix_inv;

  old_7_hecke_mtrx := HeckeOperatorDefiniteBig(M_old_7, pp);
  old_3_hecke_mtrx := HeckeOperatorDefiniteBig(M_old_3, pp);

  def_old_7_hecke_mtrxs[pp] := old_7_hecke_mtrx;
  def_old_3_hecke_mtrxs[pp] := old_3_hecke_mtrx;
end for;

primes := [pp : pp in PrimesUpTo(MAX_PRIME, F)];

for pp in primes do
  print "--------";
  print <Norm(pp), IdealOneLine(pp)>;
  print "definite";
  print Factorization(ChangeRing(CharacteristicPolynomial(def_hecke_mtrxs[pp]), F));
  print "definite new (Atkin-Lehner)";
  print Factorization(ChangeRing(CharacteristicPolynomial(def_new_hecke_mtrxs[pp]), F));

  /*
  print "definite old 7";
  print Factorization(ChangeRing(CharacteristicPolynomial(def_old_7_hecke_mtrxs[pp]), F));

  print "definite old 3";
  print Factorization(ChangeRing(CharacteristicPolynomial(def_old_3_hecke_mtrxs[pp]), F));
  */
end for;
