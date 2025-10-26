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

F := QuadraticField(2);

k := [1,2];
ideals := [N : N in IdealsUpTo(1000, F) | IsGaloisStable(N)];

for N in ideals do
  H := HeckeCharacterGroup(N, [1,2]);
  chis := [chi : chi in Elements(H) | Order(chi) eq 2 and IsCompatibleWeight(chi, k)];
  for chi in chis do
    flag := false;
    flag_quad := false;
    for psi in Elements(H) do
      if IsCompatibleWeight(psi, [1, 1]) and IsGamma1EisensteinWeight(psi, 1) then
        chi_eis := psi;
        flag := true;
        if Order(psi) eq 2 then
          flag_quad := true;
        end if;
      end if;
    end for;
    // print Norm(N), IdealOneLine(N), HeckeCharLabel(chi), chi, flag, flag_quad;
    print HeckeCharLabel(chi);
  end for;
end for;


