// for testing the degeneracy maps

load "config.m";

import "ModFrmHil/definite.m" : BasisMatrixDefinite, HeckeOperatorDefiniteBig, DegeneracyMap;
import !"Geometry/ModFrmHil/hecke.m" : restriction;
import "ModFrmHil/hecke_field.m" : DegeneracyMapDomain;

MAX_PRIME := 50;

F := QuadraticField(5);
ZF := Integers(F);

k := [4, 4];
N := 4*ZF;
pp := 2*ZF;
H := HeckeCharacterGroup(N, [1, 2]);
chi := H.0;
assert Order(chi) eq 1;

B_def := QuaternionAlgebra(1*ZF, InfinitePlaces(F) : Optimized);
O_def := MaximalOrder(B_def);

M := HilbertCuspForms(F, N, chi, k);

basis_matrix := BasisMatrixDefinite(M);
print Nrows(basis_matrix);

degen_domain := DegeneracyMapDomain(M, pp);

D1 := DegeneracyMap(degen_domain, M, 1*ZF : Big);

Up := HeckeOperatorDefiniteBig(M, pp);

D2 := D1 * Up;
