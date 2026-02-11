load "config.m";

import "ModFrmHil/definite.m" : HeckeOperatorDefiniteBig;
import !"Geometry/ModFrmHil/hecke.m" : restriction;
import !"Geometry/ModFrmHil/precompute.m" : get_rids, get_tps;

MAX_PRIME := 50;

// F := QuadraticField(6);
F := QuadraticField(29);
ZF := Integers(F);

k := [2, 2];

chi_label := "-29.0.1_91.4_2u1.0u1.2u";
chi := FullChiLabelToHeckeChar(chi_label);
N := Modulus(chi);

// N := Factorization(7*ZF)[1][1];
qq := N;
H := HeckeCharacterGroup(N, [1, 2]);
chi := H.0;
assert Order(chi) eq 1;

M := HilbertCuspForms(F, N, chi, k);

pp := PrimesUpTo(50, F)[5];
print Norm(pp);
tp := get_tps(M, pp);
print [<ki, #tp[ki]> : ki in Keys(tp)];
print [#{Norm(t) : t in tp[ki]} : ki in Keys(tp)];
A := [Rep({Norm(t) : t in tp[ki]}) : ki in Keys(tp)];


