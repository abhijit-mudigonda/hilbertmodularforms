load "config.m";

import "ModFrmHil/definite.m" : HeckeOperatorDefiniteBig;
import !"Geometry/ModFrmHil/hecke.m" : restriction;

MAX_PRIME := 50;

F := QuadraticField(2);
ZF := Integers(F);

k := [3, 3];
N := 5*ZF;

H := HeckeCharacterGroup(N, [1, 2]);
chi := H.1;
assert Order(chi) eq 4;
print HeckeCharLabel(chi);
assert HeckeCharLabel(chi) eq "-2.0.1_25.1_4u1u1.2u";

/*
B_def := QuaternionAlgebra(1*ZF, InfinitePlaces(F) : Optimized);
O_def := MaximalOrder(B_def);

M := HilbertCuspForms(F, N, chi, k);

def_hecke_mtrxs := AssociativeArray();
for pp in PrimesUpTo(MAX_PRIME, F) do
  hecke_mtrx := HeckeOperatorDefiniteBig(M, pp);
  res_hecke_mtrx := restriction(M, hecke_mtrx);
  def_hecke_mtrxs[pp] := res_hecke_mtrx;
end for;
*/

GRing := GradedRingOfHMFs(F, 400);
Mk := HMFSpace(GRing, N, k, chi);
Sk := CuspFormBasis(Mk);
Sk_hs := HeckeStabilityCuspBasis(Mk : prove:=false, stable_only:=true);
