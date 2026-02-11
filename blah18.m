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


F := QuadraticField(110);

ideals := [I : I in IdealsUpTo(2000, F) | not IsZero(I)];
for N in ideals do
  // print IdealOneLine(N);
  // print "Galois stable:", IsGaloisStable(N);
  // print #RayClassGroup(N), #RayClassGroup(N, [1]), #RayClassGroup(N, [2]), #RayClassGroup(N, [1,2]);
  if IsGaloisStable(N) and not (#RayClassGroup(N, [1,2]) gt Max([#RayClassGroup(N), #RayClassGroup(N, [1]), #RayClassGroup(N, [2])])) then
    print "!!!!!!!!!!!!! CONTRADICTION !!!!!!!!!!!!!!!!", IdealOneLine(N);
  end if;
end for;

print "number of ideals", #ideals;
