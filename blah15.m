// check the rationality condition on dihedral [1,2] forms

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
M := GradedRingOfHMFs(F, 300);
k := [1, 2];

ideals := [I : I in IdealsUpTo(3000, F) | IsGaloisStable(I)];

for N in ideals do
  print Norm(N), IdealOneLine(N);
  H := HeckeCharacterGroup(N, [1,2]);
  chis := [chi : chi in Elements(H) | IsCompatibleWeight(chi, k)];
  for chi in chis do
    if Order(chi) eq 2 then
      Mk := HMFSpace(M, N, k, chi);
      try
        if #PossibleGrossenchars(Mk) gt 0 then
          print "-----";
          print "hi!!!", chi, Order(chi), HeckeCharLabel(chi);
          print "-----";
        end if;
      catch e
        print "failed at", HeckeCharLabel(chi);
      end try;

    end if;
  end for;
e
