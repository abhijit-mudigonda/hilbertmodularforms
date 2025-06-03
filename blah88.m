load "config.m";

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

B := QuaternionAlgebra<F | -3, -12*F.1>;
O := MaximalOrder(B);
Gamma := FuchsianGroup(O);

/*
_ := Group(Gamma);

Uside := Gamma`ShimGroupSidepairs;
mside := Gamma`ShimGroupSidepairsMap;

d := #Generators(Uside);
U,m := Group(Gamma);
gquats := [m(mside(Uside.i)) : i in [1..d]];
Gamma`ShimGroupSidepairsQuats := <gquats, [g^(-1) : g in gquats]>;
*/

/*
gquats_zf := [
    [ZF.1, -1/2*ZF.1, 0, 0],
    [-ZF.2, 1/2*ZF.1 - 1/2*ZF.3, 0, -1/6*ZF.1 + 1/12*ZF.3],
    [ZF.1 - ZF.3, -ZF.1 + ZF.3, -1/3*ZF.1 + 1/6*ZF.3, 1/6*ZF.1 - 1/12*ZF.3]
];

gquats_O := [(&+[Eltseq(x)[i] * O.i : i in [1 .. 4]]) : x in gquats_zf];
gquats := [Gamma!x : x in gquats_O];


Gamma`ShimGroupSidepairsQuats := <gquats, [g^(-1) : g in gquats_O]>;
D := UnitDisc(: Center := 9/10*UpperHalfPlane().1, Precision := B`prec);
FundamentalDomain(Gamma, D : StartDomain:=gquats_O);
*/

Gamma_datum := cIdealDatum(Gamma, 3*ZF);
Gamma_datum`CosetReps;
