load "config.m";

try 
  k := [1, 2];
  chi := FullChiLabelToHeckeChar(chi_label);
  N := Modulus(chi);
  F := NumberField(Order(N));

  assert IsCompatibleWeight(chi, k);

  M := GradedRingOfHMFs(F, 440);
  Mk := HMFSpace(M, N, k, chi);
  Sk := CuspFormBasis(Mk : StableOnly:=true);

  if #Sk gt 0 then
    print "\n********************************************************************";
    print "================= Found some forms!", chi_label;
    print "********************************************************************\n";
    
    // Log labels with #Sk > 0
    nonzero_log_file := Open("logging/nonzero_sk.log", "a");
    fprintf nonzero_log_file, "%o\n", chi_label;
    delete nonzero_log_file;
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

