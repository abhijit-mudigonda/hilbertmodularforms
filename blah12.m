load "config.m";

import "ModFrmHil/definite.m" : HeckeOperatorDefiniteBig;
import !"Geometry/ModFrmHil/hecke.m" : restriction;
import !"Geometry/ModFrmHil/precompute.m" : get_rids, get_tps;

MAX_PRIME := 50;

F := QuadraticField(6);
ZF := Integers(F);

k := [2, 2];
N := Factorization(7*ZF)[1][1];
qq := N;
H := HeckeCharacterGroup(N, [1, 2]);
chi := H.0;
assert Order(chi) eq 1;


/*
Ds := [n : n in [2 .. 100] | not IsSquare(n) and (n mod 4) in [0, 1]];
for D in Ds do
  F := QuadraticField(D);
  ZF := Integers(F);
  B_def := QuaternionAlgebra(1*ZF, InfinitePlaces(F) : Optimized);
  O_def := MaximalOrder(B_def);
  print D, #RightIdealClasses(O_def);
end for;
*/

M := HilbertCuspForms(F, N, chi, k);

pp := PrimesUpTo(50, F)[5];
print Norm(pp);
tp := get_tps(M, pp);
print [<ki, #tp[ki]> : ki in Keys(tp)];
print [#{Norm(t) : t in tp[ki]} : ki in Keys(tp)];
A := [Rep({Norm(t) : t in tp[ki]}) : ki in Keys(tp)];


