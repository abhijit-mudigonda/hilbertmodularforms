load "config.m";

import "ModFrmHil/weight_rep.m" : matrix_of_induced_action;

SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);


N := 3*ZF;
k := [3,3,3];
chi := HeckeCharacterGroup(N, [1,2,3]).1;
// print HeckeCharLabel(chi);
// chi := 1;

B := QuaternionAlgebra<F | -F.1^2 + F.1 - 1, -8*F.1^2 + 4*F.1 + 16>;
O := MaximalOrder(B);
Gamma := FuchsianGroup(O);

pps := PrimesUpTo(30, F);
pp := pps[1];

X := cIdealDatum(Gamma, N : chi:=chi);
G, _, m := Group(X`FuchsianGroup);
gens := [m(G.i) : i in [1 .. #Generators(G)]];
print gens;

assert #gens eq 2;
a := matrix_of_induced_action(gens[1], k, X);
b := matrix_of_induced_action(gens[2], k, X);

IsOne(b^3);
IsOne(a^7);
IsOne((b^-1 * a^-1)^2);
