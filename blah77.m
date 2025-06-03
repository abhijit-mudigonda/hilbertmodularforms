load "config.m";

SetVerbose("ModFrmHil", 3);
SetVerbose("HilbertModularForms", 3);

F := QuadraticField(5);
ZF := Integers(F);
N := ideal<ZF | -32*ZF.2 - 19>;
k := [2,3];
H := HeckeCharacterGroup(N, [1,2]);
chi := H.2;
assert HeckeCharLabel(chi) eq "-5.0.1_55.2_2u1.0u1.2u";

M := GradedRingOfHMFs(F, 200);
Mk := HMFSpace(M, N, k, chi);
Sk := NewCuspFormBasis(Mk);

pis := [];
for pp in PrimesUpTo(50, F) do
  _, pi := IsNarrowlyPrincipal(pp);
  Append(~pis, pi);
end for;

[#Intersection(Sk, [HeckeOperatorFourier(f, F!pi) : f in Sk]) : pi in pis];


