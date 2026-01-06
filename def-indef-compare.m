load "config.m";

import "ModFrmHil/definite.m" : HeckeOperatorDefiniteBig;
import !"Geometry/ModFrmHil/hecke.m" : restriction;

MAX_PRIME := 40;


// weight [2, 2] space for sanity checking

/*
F := QuadraticField(5);
ZF := Integers(F);
k := [2, 2];
N := Factorization(7*ZF)[1][1];
qq := N;
H := HeckeCharacterGroup(N, [1, 2]);
chi := H.0;
assert Order(chi) eq 1;
*/

F := QuadraticField(2);
ZF := Integers(F);
k := [2,3];

chi_label := "-2.0.1_119.2_2u1.1u1.2u";
qq := ideal<ZF | 17, ZF.2 + 11>;

chi_label := "-2.0.1_224.1_2u1.0.1.1u1.2u";
qq := ideal<ZF | 7, ZF.2 + 3>;

chi_label := "-2.0.1_252.1_2u1.0.1u1.2u";
qq := 3*ZF;

chi_label := "-2.0.1_391.1_2u1.1u1.2u";
qq := ideal<ZF | 17, ZF.2 + 6>;

chi := FullChiLabelToHeckeChar(chi_label);
assert Order(chi) eq 2;
assert IsCompatibleWeight(chi, k);
N := Modulus(chi);
print "Norm qq", Norm(qq);
assert N ne Conductor(chi);
assert (N / qq) subset Conductor(chi);
H := HeckeCharacterGroup(N, [1,2]);

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

/*
k := [2, 2];
N := 14*ZF;
H := HeckeCharacterGroup(N, [1, 2]);
chi := H.1^2; 
assert Order(chi) eq 3;
qq := 2*ZF;
*/

// COMPUTE USING DEFINITE METHOD

print "-------- DEFINITE -----------";

B_def := QuaternionAlgebra(1*ZF, InfinitePlaces(F) : Optimized);
O_def := MaximalOrder(B_def);

// def_hm_inf := HeckeOperatorDefiniteBig(M, "Infinity");

M := HilbertCuspForms(F, N, chi, k);
M_new := NewSubspace(M, N);
_ := Dimension(M_new);  // This triggers the basis matrix computation
new_basis_matrix := M_new`basis_matrix;
new_basis_matrix_inv := M_new`basis_matrix_inv;

M_old := HilbertCuspForms(F, N / qq, Restrict(chi, N / qq, [1,2]), k);

def_hecke_mtrxs := AssociativeArray();
def_new_hecke_mtrxs := AssociativeArray();
def_old_hecke_mtrxs := AssociativeArray();
for pp in PrimesUpTo(MAX_PRIME, F) do
  hecke_mtrx := HeckeOperatorDefiniteBig(M, pp);
  // res_hecke_mtrx := restriction(M, hecke_mtrx);
  //def_hecke_mtrxs[pp] := res_hecke_mtrx;
  def_hecke_mtrxs[pp] := hecke_mtrx;
  // print <Norm(pp), IdealOneLine(pp), CharacteristicPolynomial(res_hecke_mtrx)>;
  def_new_hecke_mtrxs[pp] := new_basis_matrix * hecke_mtrx * new_basis_matrix_inv;

  old_hecke_mtrx := HeckeOperatorDefiniteBig(M_old, pp);
  def_old_hecke_mtrxs[pp] := old_hecke_mtrx;
end for;

// COMPUTE USING INDEFINITE METHOD (away from qq = 2*ZF)
/*
print "-------- INDEFINITE -----------";
B_indef := QuaternionAlgebra(qq, [InfinitePlaces(F)[1]]);
O_indef := MaximalOrder(B_indef);
Gamma := FuchsianGroup(O_indef);

_, K, _ := Splittings(B_indef);
// L := Compositum(K, CyclotomicField(Order(chi)));

// indef_hm_inf := HeckeMatrix2(Gamma, N, "Infinity", k, chi);

indef_hecke_mtrxs := AssociativeArray();

N_indef := N / qq;
chi_indef := Restrict(chi, N_indef, [1,2]);

primes := [pp : pp in PrimesUpTo(MAX_PRIME, F) | pp ne qq];
for pp in primes do
  hecke_mtrx := HeckeMatrix2(Gamma, N_indef, pp, k, chi_indef : HeckeMatrixField:=K);
  indef_hecke_mtrxs[pp] := hecke_mtrx;
  // print <Norm(pp), IdealOneLine(pp), SquareRoot(CharacteristicPolynomial(hecke_mtrx))>;
end for;

*/
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

// only use primes ne qq if you uncomment the indefinite code
primes := PrimesUpTo(MAX_PRIME, F);

for pp in primes do
  print "--------";
  print <Norm(pp), IdealOneLine(pp)>;
  print "definite";
  print Factorization(ChangeRing(CharacteristicPolynomial(def_hecke_mtrxs[pp]), F));
  // assert CharacteristicPolynomial(def_hecke_mtrxs[pp]) eq CharacteristicPolynomial(HeckeMatrix(Sk, pp));
  // assert IsDivisibleBy(CharacteristicPolynomial(HeckeMatrix(Sk, pp)), SquareRoot(CharacteristicPolynomial(indef_hecke_mtrxs[pp])));
  print "definite new";
  print Factorization(ChangeRing(CharacteristicPolynomial(def_new_hecke_mtrxs[pp]), F));

  print "definite old";
  print Factorization(ChangeRing(CharacteristicPolynomial(def_old_hecke_mtrxs[pp]), F));

//  print "indefinite";
//  print Factorization(SquareRoot(ChangeRing(CharacteristicPolynomial(indef_hecke_mtrxs[pp]), F)));
end for;

