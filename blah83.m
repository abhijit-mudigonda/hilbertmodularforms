import "ModFrmHil/diamond.m" : SaveHeckeMatrix, GetHeckeMatrix,
       SerializeAlgAssVOrdEltSeq, DeserializeAlgAssVOrdEltSeq,
       SaveXParams, UpdateIdealDatumFromSave;
load "config.m";

SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

N := ideal<ZF | [43, 0, 0], [15, 1, 0]>;
k := [2,2,3];
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.1;
assert Order(chi) eq 2;
assert HeckeCharLabel(chi) eq "1.-2.-1.1_43.2_2u1u1.2.3u";

M := HilbertCuspForms(F, N, chi, k);

O := QuaternionOrder(M);
Gamma := FuchsianGroup(O);
X := cIdealDatum(Gamma, N : chi:=chi);

X`CosetReps;
[X`ResidueMap(x) : x in X`CosetReps];
// SaveXParams(X);
/*
ord_label, gamma0cosets_ser := SerializeAlgAssVOrdEltSeq(X`CosetReps);
f := Open("blah", "w");
WriteObject(f, ord_label);
WriteObject(f, gamma0cosets_ser);

[X`ResidueMap(x) : x in X`CosetReps];
// A := X`CosetRepsByP1;
// [<ki, A[ki]> : ki in Keys(A)];

// SaveXandk(Gamma, N, chi, k);
*/
