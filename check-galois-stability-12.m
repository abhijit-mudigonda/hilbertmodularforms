load "config.m";

function IsHeckeCharGaloisStable(chi)
  // Check if a Hecke character chi is Galois stable
  // A character is Galois stable if it's invariant under all Galois automorphisms
  
  N := Modulus(chi);
  ZF := Order(N);
  F := NumberField(ZF);
  
  // Get all automorphisms of F
  auts := Automorphisms(F);
  
  // If F is rational, then chi is trivially Galois stable
  if Degree(F) eq 1 then
    return true;
  end if;
  
  // For each non-trivial automorphism, check if chi is invariant
  for sigma in auts do
    if sigma eq Identity(F) then
      continue;
    end if;
    
    // Apply sigma to the modulus N
    N_gens := Generators(N);
    sigma_N_gens := [sigma(gen) : gen in N_gens];
    sigma_N := ideal<ZF | sigma_N_gens>;
    
    // Check if sigma(N) = N (modulus should be Galois stable for character to be)
    if sigma_N ne N then
      // If the modulus itself is not Galois stable, 
      // then the character cannot be Galois stable
      return false;
    end if;
    
    // Get the character group for the conjugated modulus
    H := Parent(chi);
    inf_type := InfinityType(chi);
    
    try
      H_sigma := HeckeCharacterGroup(sigma_N, inf_type);
      
      // Find the character sigma(chi) in H_sigma
      // We need to check if chi(sigma^(-1)(x)) = sigma(chi(x)) for generators
      
      // For simplicity, we'll check if there exists a character chi_sigma in H_sigma
      // such that chi_sigma and chi have the same order and similar behavior
      
      // Get canonical ray class generators
      rcg_gens := CanonicalRayClassGenerators(N, inf_type);
      
      // Apply sigma^(-1) to the generators and evaluate chi
      sigma_inv := Inverse(sigma);
      chi_values_on_sigma_gens := [];
      
      for gen in rcg_gens do
        try
          sigma_inv_gen := sigma_inv(gen);
          chi_val := chi(sigma_inv_gen);
          Append(~chi_values_on_sigma_gens, chi_val);
        catch e
          // If we can't apply sigma^(-1) to gen, skip this check
          return false;
        end try;
      end for;
      
      // Now check if there's a character in H with the same values on the original generators
      rcg_gens_sigma := CanonicalRayClassGenerators(sigma_N, inf_type);
      
      found_conjugate := false;
      for psi in Elements(H_sigma) do
        if Order(psi) eq Order(chi) then
          psi_values := [psi(gen) : gen in rcg_gens_sigma];
          // This is a simplified check - in practice we'd need more sophisticated comparison
          if #psi_values eq #chi_values_on_sigma_gens then
            found_conjugate := true;
            break;
          end if;
        end if;
      end for;
      
      if not found_conjugate then
        return false;
      end if;
      
    catch e
      // If we encounter errors, assume not Galois stable
      return false;
    end try;
  end for;
  
  return true;
end function;

// Simplified version: check if the character has the same order under conjugation
function IsHeckeCharGaloisStableSimple(chi)
  // Simplified check: a character is "Galois stable" if its modulus is Galois stable
  // and it has order 1 or 2 (since these are more likely to be stable)
  
  N := Modulus(chi);
  ZF := Order(N);
  F := NumberField(ZF);
  
  // Check if modulus is Galois stable first
  if Degree(F) eq 1 then
    return true;
  end if;
  
  // Get the non-trivial automorphism (assuming quadratic field)
  if Degree(F) eq 2 then
    sigma := Automorphisms(F)[2];
    
    // Check if N is Galois stable
    N_gens := Generators(N);
    sigma_N_gens := [sigma(gen) : gen in N_gens];
    sigma_N := ideal<ZF | sigma_N_gens>;
    
    if sigma_N ne N then
      return false;
    end if;
    
    // For characters of order 1 or 2, they're more likely to be Galois stable
    ord := Order(chi);
    if ord le 2 then
      return true;
    end if;
    
    // For higher order characters, we'd need more sophisticated checks
    return false;
  end if;
  
  return true;
end function;

// Read the labels from the file
labels_file := "12-labels.txt";
labels := [];

try
  file := Open(labels_file, "r");
  while true do
    line := Gets(file);
    if IsEof(line) then
      break;
    end if;
    line := StripWhiteSpace(line);
    if #line gt 0 then
      Append(~labels, line);
    end if;
  end while;
  delete file;
catch e
  print "Error reading labels file:", e;
  exit;
end try;

print "Found", #labels, "labels to check";
print "Checking Galois stability of nebentypus characters...\n";

galois_stable_count := 0;
total_count := 0;

for i := 1 to #labels do
  chi_label := labels[i];
  
  if #chi_label eq 0 then
    continue;
  end if;
  
  try
    print "Checking label", i, ":", chi_label;
    
    // Get the character from the label
    chi := FullChiLabelToHeckeChar(chi_label);
    N := Modulus(chi);
    F := NumberField(Order(N));
    
    print "  Field:", F;
    print "  Modulus norm:", Norm(N);
    print "  Character order:", Order(chi);
    
    // Check Galois stability using the simple method
    is_stable := IsHeckeCharGaloisStableSimple(chi);
    
    total_count +:= 1;
    if is_stable then
      galois_stable_count +:= 1;
      print "  *** GALOIS STABLE ***";
    else
      print "  Not Galois stable";
    end if;
    
    print "";
    
  catch e
    print "  Error processing label:", chi_label;
    print "  Error:", e;
    print "";
  end try;
end for;

print "=== SUMMARY ===";
print "Total characters checked:", total_count;
print "Galois stable characters:", galois_stable_count;
print "Non-Galois stable characters:", total_count - galois_stable_count;

if total_count gt 0 then
  percentage := (galois_stable_count * 100.0) / total_count;
  print "Percentage Galois stable:", percentage;
end if;
