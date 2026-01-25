// this is code to try to actually compute what the field E (the one over which
// the tensor induction is defined) in some example. I don't remember now but I think
// it didn't work. I guess this is indication that this (CM) form isn't actually
// an example of what we want? I don't remember, should check everything carefully.

load "config.m";

// SetVerbose("ModFrmHil", 3);
PREC := 777;

// chi_label := "-5.0.1_576.1_2u0.1.0.1u1.2u";
// chi_label := "-5.0.1_784.1_2u0.1.0u1.2u";
// chi_label := "-5.0.1_1024.1_2u1.1.1u1.2u";
chi_label := "-5.0.1_1280.1_2u0.0.1.1u1.2u";

chi := FullChiLabelToHeckeChar(chi_label);
assert Order(chi) eq 2;
N := Modulus(chi);
F := NumberField(Order(N));
ZF := Integers(F);
M := GradedRingOfHMFs(F, PREC);
k := [1,2];
Mk := HMFSpace(M, N, k, chi);
Sk := CuspFormBasis(Mk);
assert #Sk eq 2;

/*
M24 := HMFSpace(M, N, [2,4]);
S24 := CuspFormBasis(M24);
assert #LinearDependence(S24) eq 0;

Sk_squared := [Sk[1]^2, Sk[1] * Sk[2], Sk[2]^2];
assert #Intersection(Sk_squared, S24) eq 3;
*/

eigs := Eigenbasis(Mk, Sk : P:=12, coprime_only:=false);
dinv_F := IdealToRep(M, 1*ZF);
aut := Automorphisms(F)[2];

for f in eigs do
  K := CoefficientRing(f);

  f_div := f / Coefficient(f, 1*ZF, dinv_F);

  /*
  for pp in PrimesUpTo(50, F) do
    _, pi := IsNarrowlyPrincipal(pp);
    assert (Coefficient(f, 1*ZF, pi * dinv_F) / Coefficient(f, 1*ZF, dinv_F)) eq Coefficient(f_div, 1*ZF, pi * dinv_F);
  end for;
  */

  flag := true;
  for pp in PrimesUpTo(100, F) do
    if IsPrime(Norm(pp)) and Norm(pp) ne 2 then 
      _, pi := IsNarrowlyPrincipal(pp);
      a_pi := Coefficient(f_div, 1*ZF, pi * dinv_F);
      a_autpi := Coefficient(f_div, 1*ZF, aut(pi) * dinv_F);
      a := a_pi * a_autpi;
      if not IsCoercible(Rationals(), a) then
        flag := false;
      end if;
      print Norm(pp), IsCoercible(Rationals(), a), MinimalPolynomial(a);
    end if;
  end for;
  print "-------------";
end for;

/*
[#LinearDependence([f, HeckeOperator(f, pp)]) : pp in PrimesUpTo(20, F)];
#Intersection([f^2], S24);
*/
