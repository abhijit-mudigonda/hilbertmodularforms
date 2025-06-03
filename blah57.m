load "config.m";

F := QuadraticField(5);
ZF := Integers(F);
M := GradedRingOfHMFs(F, 500);
N := Factorization(11*ZF)[1][1];
k := [2,4];

Mk := HMFSpace(M, N, k);
Sk := CuspFormBasis(Mk);

f := Sk[1];
f := DivideByFirstNonzeroIdlCoeff(f);
nus := Keys(M`RepToIdeal[1*ZF]);
nus := SetToSequence(nus);

[IsIntegral(Coefficient(f, 1*ZF, nu)) : nu in nus];



