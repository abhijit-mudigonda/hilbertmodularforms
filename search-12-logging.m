load "config.m";

F := QuadraticField(2);
ZF := Integers(F);
M := GradedRingOfHMFs(F, 777);
k := [1, 2];

ideals := [N : N in IdealsUpTo(600, F) | Norm(N) eq StringToInteger(N_norm)];

for N in ideals do
  H := HeckeCharacterGroup(N, [1,2]);
  for chi in Elements(H) do
    if IsCompatibleWeight(chi, k) and Order(chi) eq 2 then
      chi_label := HeckeCharLabel(chi);
      print "Trying:", IdealOneLine(N), chi_label, chi;
      cond_fin, cond_inf := Conductor(chi);
      Mk := HMFSpace(M, N, k, chi);
      try
        Sk := CuspFormBasis(Mk : StableOnly:=true);
        if #Sk gt 0 then
          print "******************** Found a [1,2] form:", IdealOneLine(N), chi, chi_label;
          // Log successful cases to found_12.log
          log_file := Open("found_12.log", "a");
          fprintf log_file, "%o %o %o\n", IdealOneLine(N), chi, chi_label;
          delete log_file;
        end if;
      catch e;
        // Extract error message as string
        error_str := Sprint(e);
        
        // Check if this is the specific Eisenstein error
        if "No Eisenstein series seems to work?" in error_str then
          // Log to no_eisenstein.log
          log_file := Open("no_eisenstein.log", "a");
          fprintf log_file, "%o %o %o\n", IdealOneLine(N), chi, chi_label;
          delete log_file;
        else
          // Log to errors.log for other errors
          log_file := Open("errors.log", "a");
          fprintf log_file, "%o %o %o\n", IdealOneLine(N), chi, chi_label;
          delete log_file;
        end if;
      end try;
    end if;
  end for;
end for;
