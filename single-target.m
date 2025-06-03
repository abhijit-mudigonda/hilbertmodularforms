import "ModFrmHil/diamond.m" : SaveHeckeMatrix; 
load "config.m";
// SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

B := QuaternionAlgebra<F | -3, -12*F.1>;
O := MaximalOrder(B);
Gamma := FuchsianGroup(O);

// _ := cIdealDatum(Gamma, N);

N := ideal<ZF | [43, 0, 0], [15, 1, 0]>;
k := [2,2,3];
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.1;
assert Order(chi) eq 2;
assert HeckeCharLabel(chi) eq "1.-2.-1.1_43.2_2u1u1.2.3u";

/*
N := 3*ZF;
k := [2,2,4];
chi := HeckeCharacterGroup(N, [1,2,3]).0;
*/

/*
N := 3*ZF;
k := [2,2,2];
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.0;
*/

/*
N := 12*ZF;
k := [2,2,3];
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.4;
assert Order(chi) eq 2;
assert HeckeCharLabel(chi) eq "1.-2.-1.1_1728.1_2u1.1.1.1u1.2.3u";
*/

CORES := StringToInteger(CORES);
CORE := StringToInteger(CORE);
MIN_NORM := StringToInteger(MIN_NORM);
MAX_NORM := StringToInteger(MAX_NORM);

pps := [pp : pp in PrimesUpTo(MAX_NORM, F) | Norm(pp) ge MIN_NORM];
pps := [pps[1 + CORES*i + CORE] : i in [0 .. Floor((#pps - 1 - CORE) / CORES)]];
print "pps", [Norm(pp) : pp in pps];

if CORE eq 0 then
  SaveHeckeMatrix(Gamma, N, "Infinity", k, chi);
end if;

for pp in pps do
  SaveHeckeMatrix(Gamma, N, pp, k, chi);
end for;
