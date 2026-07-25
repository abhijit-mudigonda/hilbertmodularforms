load "config.m";

PREC := 1000;
F := QuadraticField(5);
M := GradedRingOfHMFs(F, PREC);

try
  N := LMFDBIdeal(F, N_label);
  Mk := HMFSpace(M, N, [2,2]);
  _ := NewCuspFormBasis(Mk : SaveAndLoad:=true);
catch e
  log_file := Open("logging/precompute_errors.log", "a");
  fprintf log_file, "%o %o\n", N_label, Sprint(e);
  delete log_file;
end try;

quit;
