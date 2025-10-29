load "config.m";

SetVerbose("ModFrmHil", 3);
PREC := 500;
F := QuadraticField(2);
ZF := Integers(F);

M := GradedRingOfHMFs(F, PREC);

full_label := "-2.0.1_119.2_2u1.0u1.2u";

chi := FullChiLabelToHeckeChar(full_label);
assert Order(chi) eq 2;
N := Modulus(chi);
k := [1,2];
Mk := HMFSpace(M, N, k, chi);
Sk := CuspFormBasis(Mk);
assert #Sk eq 2;

M24 := HMFSpace(M, N, [2,4]);
S24 := CuspFormBasis(M24);
assert #LinearDependence(S24) eq 0;

Sk_squared := [Sk[1]^2, Sk[1] * Sk[2], Sk[2]^2];
assert #Intersection(Sk_squared, S24) eq 3;

eigs := Eigenbasis(Mk, Sk : P:=12, coprime_only:=false);
dinv_F := IdealToRep(M, 1*ZF);
f := eigs[1];
K := CoefficientRing(f);
f := f / Coefficient(f, 1*ZF, dinv_F);


[#LinearDependence([f, HeckeOperator(f, pp)]) : pp in PrimesUpTo(20, F)];

#Intersection([f^2], S24);

