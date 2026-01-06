// for testing the degeneracy maps

load "config.m";

import "ModFrmHil/definite.m" : BasisMatrixDefinite, HeckeOperatorDefiniteBig, DegeneracyMap, AtkinLehnerDefiniteBig;
import !"Geometry/ModFrmHil/hecke.m" : restriction;
import "ModFrmHil/hecke_field.m" : DegeneracyMapDomain;

MAX_PRIME := 50;

F := QuadraticField(5);
ZF := Integers(F);

k := [2, 2];
N := 11*ZF;
H := HeckeCharacterGroup(N, [1, 2]);
chi := H.0;
assert Order(chi) eq 1;
pp := 2*ZF;

B_def := QuaternionAlgebra(1*ZF, InfinitePlaces(F) : Optimized);
O_def := MaximalOrder(B_def);

M_old := HilbertCuspForms(F, N, chi, k);
M := HilbertCuspForms(F, N * pp, chi, k);

basis_matrix_old := BasisMatrixDefinite(M_old);
print "level N basis dim", Nrows(basis_matrix_old);

basis_matrix := BasisMatrixDefinite(M);
print "level N * pp basis dim", Nrows(basis_matrix);

degen_domain := DegeneracyMapDomain(M, Level(M) / pp);
D1 := DegeneracyMap(degen_domain, M, 1*ZF : Big, EisensteinAllowed:=false);
D2 := DegeneracyMap(degen_domain, M, pp : Big, EisensteinAllowed:=false);
AL_2 := AtkinLehnerDefiniteBig(M, pp);
q1 := Factorization(11*ZF)[1][1];
q2 := Factorization(11*ZF)[2][1];
AL_q1 := AtkinLehnerDefiniteBig(M, q1);
AL_q2 := AtkinLehnerDefiniteBig(M, q2);

Up := HeckeOperatorDefiniteBig(M, pp);

old_hms := AssociativeArray();
hms := AssociativeArray();

MAX_PRIME := 70;
for qq in PrimesUpTo(MAX_PRIME, F) do
  old_hms[qq] := HeckeOperatorDefiniteBig(M_old, qq);
  hms[qq] := HeckeOperatorDefiniteBig(M, qq);
end for;

v := KernelMatrix(old_hms[2*ZF]);

Dold := VerticalJoin(D1, D2);

AL_q1_res := Solution(D1, D1 * AL_q1);
AL_q2_res := Solution(D1, D1 * AL_q2);

AL_2_res := Solution(Dold, Dold * AL_2);

for qq in PrimesUpTo(MAX_PRIME, F : coprime_to:=22) do
  hm_res := Solution(D1, D1 * hms[qq]);
  print Norm(qq), Eigenvalues(hm_res);
end for;

print "---";
for qq in PrimesUpTo(MAX_PRIME, F : coprime_to:=22) do
  hm_res := Solution(D2, D2 * hms[qq]);
  print Norm(qq), Eigenvalues(hm_res);
end for;

assert Nrows(D1) eq 3;

// Check the relation Up ∘ β = α in the non-Big representation
print "\n=== Checking Up ∘ β = α ===";

// Compute non-Big degeneracy maps
// α: M_old -> M (the degeneracy map by 1)
// β: M_old -> M (the degeneracy map by pp)
alpha := DegeneracyMap(degen_domain, M, 1*ZF : Big:=false, EisensteinAllowed:=false);
beta := DegeneracyMap(degen_domain, M, pp : Big:=false, EisensteinAllowed:=false);

print "alpha dimensions:", Nrows(alpha), "x", Ncols(alpha);
print "beta dimensions:", Nrows(beta), "x", Ncols(beta);

// Compute Up restricted to M (in non-Big representation)
// This is the Hecke operator U_pp acting on M
Up_restricted := Solution(basis_matrix, basis_matrix * Up);
print "Up_restricted dimensions:", Nrows(Up_restricted), "x", Ncols(Up_restricted);

// Check if β * Up_restricted = α
// This should hold if Up ∘ β = α
print "hoi", beta * Up_restricted eq alpha;
print "hoii", alpha * Up_restricted eq beta;

