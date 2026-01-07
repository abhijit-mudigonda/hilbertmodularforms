// this is code to try to actually compute what the field E (the one over which
// the tensor induction is defined) in some example. I don't remember now but I think
// it didn't work. I guess this is indication that this (CM) form isn't actually
// an example of what we want? I don't remember, should check everything carefully.

load "config.m";

SetVerbose("ModFrmHil", 3);
PREC := 777;
// F := QuadraticField(2);
F = QuadraticField(5);
ZF := Integers(F);

M := GradedRingOfHMFs(F, PREC);

// full_label := "-2.0.1_119.2_2u1.0u1.2u";
full_label := "-5.0.1_576.1_2u0.1.0.1u1.2u";

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

dinv_F := IdealToRep(M, 1*ZF);

f_div := f / Coefficient(f, 1*ZF, dinv_F);
aut := Automorphisms(F)[2];
for pp in PrimesUpTo(50, F) do
  _, pi := IsNarrowlyPrincipal(pp);
  assert (Coefficient(f, 1*ZF, pi * dinv_F) / Coefficient(f, 1*ZF, dinv_F)) eq Coefficient(f_div, 1*ZF, pi * dinv_F);
end for;

for pp in PrimesUpTo(100, F) do
  if IsPrime(Norm(pp)) and Norm(pp) ne 2 then 
    _, pi := IsNarrowlyPrincipal(pp);
    a_pi := Coefficient(f_div, 1*ZF, pi * dinv_F);
    a_autpi := Coefficient(f_div, 1*ZF, aut(pi) * dinv_F);
    a := a_pi * a_autpi;
    print Norm(pp), IsCoercible(Rationals(), a), Degree(MinimalPolynomial(a));
  end if;
end for;

[#LinearDependence([f, HeckeOperator(f, pp)]) : pp in PrimesUpTo(20, F)];

#Intersection([f^2], S24);
