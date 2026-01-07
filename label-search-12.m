load "config.m";

square_check := true;

try 
  k := [1, 2];
  chi := FullChiLabelToHeckeChar(chi_label);
  N := Modulus(chi);
  F := NumberField(Order(N));

  assert IsCompatibleWeight(chi, k);

  M := GradedRingOfHMFs(F, 300);
  Mk := HMFSpace(M, N, k, chi);
  Sk := CuspFormBasis(Mk : StableOnly:=true);

  pp_2 := PrimesUpTo(20, F)[2];
  if #Sk gt 0 then
    V := HeckeStableSubspace(Sk, pp_2 : use_fourier:=true);
  end if;

  if #Sk gt 0 and #V gt 0 then
    print "Found some forms, need to verify", chi_label;
        
    Sk := V;
    // Log labels with #Sk > 0
    nonzero_log_file := Open("logging/nonzero_sk.log", "a");
    fprintf nonzero_log_file, "%o\n", chi_label;
    delete nonzero_log_file;

    Sk_squared := &cat[[Sk[i] * Sk[j] : j in [1 .. #Sk] | j ge i] : i in [1 .. #Sk]];
    
    M24 := HMFSpace(M, N, [2,4]);
    S24 := CuspFormBasis(M24);
    if #LinearDependence(S24) eq 0 then
      if #LinearDependence(Sk_squared) gt 0 then
        Sk_squared := Basis(Sk_squared);
      end if;
      if #Intersection(Sk_squared, S24) eq #Sk_squared then
        print "\n********************************************************************";
        print "================= Found some forms!", chi_label;
        print "********************************************************************\n";
        proven_nonzero_log_file := Open("logging/squarechecked_nonzero_sk.log", "a");
        fprintf proven_nonzero_log_file, "%o\n", chi_label;
        delete proven_nonzero_log_file;
      end if;
    end if;
  else
    print "No forms at", chi_label;
    // Log labels with #Sk = 0
    zero_log_file := Open("logging/zero_sk.log", "a");
    fprintf zero_log_file, "%o\n", chi_label;
    delete zero_log_file;
  end if;
catch e
  // Log errors
  error_str := Sprint(e);
  print "Error occurred:", error_str;
  
  log_file := Open("logging/errors.log", "a");
  fprintf log_file, "%o %o", chi_label, error_str;
  delete log_file;
end try;
