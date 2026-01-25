load "config.m";

function IsGaloisStable(N)
  // Check if an integral ideal N of a quadratic field is Galois stable
  ZF := Order(N);
  F := NumberField(ZF);
  
  // Get the non-trivial automorphism of F
  sigma := Automorphisms(F)[2];
  
  // Factor the ideal N
  factorization := Factorization(N);
  
  // Check each prime factor
  for factor in factorization do
    p := factor[1];  // the prime ideal
    e := factor[2];  // the exponent
    
    // Compute the conjugate of the prime ideal p under sigma
    p_gens := Generators(p);
    sigma_p_gens := [sigma(gen) : gen in p_gens];
    sigma_p := ideal<ZF | sigma_p_gens>;
    
    // Check if p is split (i.e., p != sigma(p))
    if p ne sigma_p then
      // p is split, so check if sigma(p) appears with the same exponent
      if Valuation(N, sigma_p) ne e then
        return false;
      end if;
    end if;
  end for;
  
  return true;
end function;


F := QuadraticField(5);
ZF := Integers(F);
k := [2,3];
M := GradedRingOfHMFs(F, 444);
dinv_F := IdealToRep(M, 1*ZF);
aut := Automorphisms(F)[2];

for N in [I : I in IdealsUpTo(250, F) | Norm(I) gt 140] do
  H := HeckeCharacterGroup(N, [1,2]);
  chis := [chi : chi in Elements(H) | IsCompatibleWeight(chi, k) and Order(chi) eq 2];
  for chi in chis do
    Mk := HMFSpace(M, N, k, chi);
    Sk := CuspFormBasis(Mk);
    if #Sk gt 0 then
      print "Here be forms", IdealOneLine(N), chi;
      eigs := Eigenbasis(Mk, Sk : P:=12);
      for f in eigs do
        K := CoefficientRing(f);

        a_1 := Coefficient(f, 1*ZF, dinv_F);
        if not IsZero(a_1) then
          f_div := f / a_1;
          flag := true;
          for pp in PrimesUpTo(300, F) do
            if IsPrime(Norm(pp)) and Norm(pp) ne 2 then 
              _, pi := IsNarrowlyPrincipal(pp);
              a_pi := Coefficient(f_div, 1*ZF, pi * dinv_F);
              a_autpi := Coefficient(f_div, 1*ZF, aut(pi) * dinv_F);
              a := a_pi * a_autpi;
              if not IsCoercible(Rationals(), a) then
                flag := false;
              end if;
              print Norm(pp), IsCoercible(Rationals(), a), Degree(MinimalPolynomial(a));
            elif not IsPrime(Norm(pp)) then
              _, pi := IsNarrowlyPrincipal(pp);
              a_pi := Coefficient(f_div, 1*ZF, pi * dinv_F);
              print Norm(pp), IsCoercible(Rationals(), a_pi), Degree(MinimalPolynomial(a_pi));
            end if;
          end for;
          if flag then
            print "Found one that works", IdealOneLine(N), chi, HeckeCharLabel(chi);
          end if;
        end if;
        print "--------";
      end for;
    end if;
  end for;
end for;

