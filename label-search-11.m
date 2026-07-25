load "config.m";

PREC := assigned prec select StringToInteger(prec) else 1000;
CLASSIFY_PRIME_BOUND := 100;     // primes (coprime to N) checked for order-4/5 Frobenius elements
EIGENBASIS_PRIME_BOUND := 20;    // matches the P:=20 used in a4-form.m's Eigenbasis call
TOL := 10^-6;                    // numeric tolerance for the orderer

function orderer(x)
  targets := [<RealField(30)!4, 1>, <RealField(30)!0, 2>, <RealField(30)!1, 3>,
              <RealField(30)!2, 4>, <RealField(30)!3, 6>,
              <(3 + Sqrt(RealField(30)!5))/2, 5>, <(3 - Sqrt(RealField(30)!5))/2, 5>];
  for c in Conjugates(x) do
    if Abs(Im(c)) lt TOL then
      re := Re(c);
      for t in targets do
        if Abs(re - t[1]) lt TOL then
          return t[2];
        end if;
      end for;
    end if;
  end for;
  return -1; // unrecognized - log for manual inspection instead of crashing the batch
end function;

try
  chi := FullChiLabelToHeckeChar(chi_label);
  N := Modulus(chi);
  F := NumberField(Order(N));
  k := [1, 1];
  assert IsCompatibleWeight(chi, k);

  M := GradedRingOfHMFs(F, PREC);
  Mk := HMFSpace(M, N, k, chi);
  // Calling HeckeStabilityCuspBasis directly (rather than
  // CuspFormBasis(Mk : StableOnly:=true, SaveAndLoad:=true)) is required to get
  // both stable_only and SaveAndLoad honored at once: CuspFormBasis's top-level
  // SaveAndLoad branch ignores StableOnly and always does the full
  // (slower) holomorphy-proving pass. Calling this directly also means
  // SaveAndLoad only caches the weight-2 "upstairs" space and its divisors
  // (where the actual repeated work across labels sharing a level lives),
  // not the weight-1 result itself.
  Sk := HeckeStabilityCuspBasis(Mk : prove:=false, stable_only:=true, SaveAndLoad:=true);
  Ek := EisensteinBasis(Mk);

  if #Sk - #Intersection(Sk, Ek) gt 0 then
    log_file := Open("logging/11_nonzero.log", "a");
    fprintf log_file, "%o\n", chi_label;
    delete log_file;

    B := (#Sk gt 1) select Eigenbasis(Mk, Sk : P:=EIGENBASIS_PRIME_BOUND) else Sk;
    B := [DivideByFirstNonzeroIdlCoeff(f) : f in B];

    classify_primes := PrimesUpTo(CLASSIFY_PRIME_BOUND, F : coprime_to:=N);
    for f in B do
      if #Intersection([f], Ek) eq 0 and not ProbabilisticDihedralTest(f) then
        orders := [orderer(Coefficient(f, pp)^2 / chi(pp)) : pp in classify_primes];

        log_name := "logging/11_a4_candidate.log";
        if 4 in orders then
          log_name := "logging/11_s4_found.log";
        elif 5 in orders then
          log_name := "logging/11_a5_found.log";
        elif -1 in orders then
          log_name := "logging/11_unrecognized.log";
        end if;

        log_file := Open(log_name, "a");
        fprintf log_file, "%o %o\n", chi_label, orders;
        delete log_file;
      end if;
    end for;
  end if;
catch e
  error_str := Sprint(e);
  log_name := ("No Eisenstein series seems to work?" in error_str)
    select "logging/11_no_eisenstein.log" else "logging/11_errors.log";
  log_file := Open(log_name, "a");
  fprintf log_file, "%o %o\n", chi_label, error_str;
  delete log_file;
end try;

quit;
