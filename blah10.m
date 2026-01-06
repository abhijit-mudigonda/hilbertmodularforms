// for testing the degeneracy maps

load "config.m";

import "ModFrmHil/definite.m" : BasisMatrixDefinite, HeckeOperatorDefiniteBig, DegeneracyMap, AtkinLehnerDefiniteBig;
import !"Geometry/ModFrmHil/hecke.m" : restriction;
import "ModFrmHil/hecke_field.m" : DegeneracyMapDomain;

MAX_PRIME := 50;

F := QuadraticField(5);
ZF := Integers(F);

k := [4, 4];
N := 6*ZF;
pp := 3*ZF;
H := HeckeCharacterGroup(N, [1, 2]);
chi := H.0;
assert Order(chi) eq 1;

/*
F := QuadraticField(2);
ZF := Integers(F);
N := ideal<ZF | 119, -ZF.2 - 45>;
k := [2, 3];
H := HeckeCharacterGroup(N, [1, 2]);
chi := H.2;
assert IsCompatibleWeight(chi, [2,3]);
pp := Conductor(chi);
assert Norm(pp) eq 7;
*/

B_def := QuaternionAlgebra(1*ZF, InfinitePlaces(F) : Optimized);
O_def := MaximalOrder(B_def);

M := HilbertCuspForms(F, N, chi, k);
M_old := HilbertCuspForms(F, N / pp, Restrict(chi, N / pp, [1 .. Degree(F)]), k);
M_new := HilbertCuspForms(F, N * pp, Extend(chi, N * pp, [1 .. Degree(F)]), k);

basis_matrix_old := BasisMatrixDefinite(M_old);
print "level N/pp basis dim", Nrows(basis_matrix_old);

basis_matrix := BasisMatrixDefinite(M);
print "level N basis dim", Nrows(basis_matrix);

basis_matrix_new := BasisMatrixDefinite(M_new);
print "level N * pp basis dim", Nrows(basis_matrix_new);

degen_domain := DegeneracyMapDomain(M, Level(M) / pp);
D1 := DegeneracyMap(degen_domain, M, 1*ZF : Big);
D2 := DegeneracyMap(degen_domain, M, pp : Big);
AL := AtkinLehnerDefiniteBig(M, pp);
Up := HeckeOperatorDefiniteBig(M, pp);

AL_res := restriction(M, AL);
Up_res := restriction(M, Up);

// prints true
print Rowspace(D1 * Up) eq Rowspace(D2);
print Rowspace(D2 * Up) eq Rowspace(D1);

// prints false!
print Rowspace(D1 * AL) eq Rowspace(D2);

degen_domain_new := DegeneracyMapDomain(M_new, Level(M_new) / pp);
D1_new := DegeneracyMap(degen_domain_new, M_new, 1*ZF : Big);
D2_new := DegeneracyMap(degen_domain_new, M_new, pp : Big);
Up_new := HeckeOperatorDefiniteBig(M_new, pp);
AL_new := AtkinLehnerDefiniteBig(M_new, pp);

AL_new_res := restriction(M_new, AL_new);
Up_new_res := restriction(M_new, Up_new);

print Rowspace(D1_new * Up_new) eq Rowspace(D2_new);
print Rowspace(D2_new * Up_new) eq Rowspace(D1_new);

// prints false!
print Rowspace(D1_new * AL_new) eq Rowspace(D2_new);
