import "ModFrmHil/diamond.m" : SaveHeckeMatrix, GetHeckeMatrix,
       SaveXandk, XandkFilepath, SerializeAlgAssVOrdEltSeq;
load "config.m";

SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

/*
N := ideal<ZF | [43, 0, 0], [15, 1, 0]>;
k := [2,2,3];
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.1;
assert Order(chi) eq 2;
assert HeckeCharLabel(chi) eq "1.-2.-1.1_43.2_2u1u1.2.3u";
*/

N := 3*ZF;
k := [2,2,2];

M := HilbertCuspForms(F, N, chi, k);

O := QuaternionOrder(M);
Gamma := FuchsianGroup(O);
X := cIdealDatum(Gamma, N : chi:=chi);

ord_label, gamma0cosets_ser := SerializeAlgAssVOrdEltSeq(X`CosetReps);

savefile_name := XFilepath(Gamma, N, chi);
savefile := Open(savefile_name, "w");
WriteObject(savefile, ord_label);
WriteObject(savefile, gamma0cosets_ser);
WriteObject(savefile, gamma
savefile := 0;

