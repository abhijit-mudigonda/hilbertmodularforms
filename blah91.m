load "config.m";

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

B := QuaternionAlgebra<F | -3, -12*F.1>;
O := MaximalOrder(B);
Gamma := FuchsianGroup(O);

N := 3*ZF;
// _ := cIdealDatum(Gamma, N);

k := [2,2,2];
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.0;

// THE FOLLOWING BLOCKS GET INSERTED HERE

save "blah91-save"; 
restore "blah91-save";

// For some stupid reason, there will be a seg fault
// if some commands are executed immediately after a restore.
// This includes the Sleep() command, which I'd usually use to
// wait. Anyways, this is the only way I could find to kill time
// that doesn't also cause a seg fault.

for i in [1 .. 100000000] do
  j := i;
end for;

pp := Factorization(7*ZF)[1][1];
HeckeMatrix2(Gamma, N, pp, k, chi);
