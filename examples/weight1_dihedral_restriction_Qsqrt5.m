load "config.m";

/***************************************************************
* Studies the restriction to the diagonal of parallel weight
* [1,1] dihedral Hilbert modular eigenforms over F = Q(sqrt(5)),
* across levels N of increasing norm.
*
* For each level N and each nebentypus chi compatible with
* weight [1,1], we compute the dihedral subspace DihedralBasis(Mk)
* (whose elements are already Hecke eigenforms -- see ThetaSeries
* in ModFrmHilD/Creation/Dihedral.m), restrict each eigenform to
* the diagonal, and record the resulting classical q-expansion
* coefficients. We use RestrictionToDiagonal(... : AsCoefficients
* := true) so that we get the raw coefficients directly (in
* whatever ring they live in) instead of a classical ModFrmElt,
* since the latter path only supports rational coefficients.
*
* Output:
*   WeightOneDihedralRestriction/summary.csv       -- per-level counts
*   WeightOneDihedralRestriction/qexpansions/*.txt  -- one file per
*     (level, character, eigenform, component) with the q-expansion
*     coefficients of the restriction
*   WeightOneDihedralRestriction/errors.log         -- any failures
*     encountered while building a dihedral subspace
***************************************************************/

F := QuadraticField(5);
ZF := Integers(F);
k := [1, 1];
n := Degree(F);

MAX_LEVEL_NORM := 600;   // norm(N) search bound -- first dihedral form is at norm 241
HMF_PRECISION  := 500;   // Fourier precision of the ambient graded ring; empirically
                         // yields ~21 restriction coefficients (coefficient count
                         // grows like sqrt(HMF_PRECISION), see coeff_bound_scan.m)
NUM_COEFFS     := 20;    // number of classical q-expansion coefficients to save

OUT_DIR := "WeightOneDihedralRestriction";
QEXP_DIR := OUT_DIR * "/qexpansions";
System("mkdir -p " * QEXP_DIR);

summary_fp := Open(OUT_DIR * "/summary.csv", "w");
Puts(summary_fp, "level_norm,level_label,level_oneline,num_characters,num_dihedral_eigenforms,num_zero_restriction");
Flush(summary_fp);

errors_fp := Open(OUT_DIR * "/errors.log", "w");

M := GradedRingOfHMFs(F, HMF_PRECISION);

ideals := [N : N in IdealsUpTo(MAX_LEVEL_NORM, F) | not IsZero(N)];
printf "Studying weight [1,1] dihedral restrictions for %o levels up to norm %o.\n", #ideals, MAX_LEVEL_NORM;

for N in ideals do
  level_label := Join(Split(LMFDBLabel(N), "."), "_");
  H := HeckeCharacterGroup(N, [1,2]);
  chis := [chi : chi in Elements(H) | IsCompatibleWeight(chi, k)];

  num_dihedral := 0;
  num_zero := 0;

  for chi_idx -> chi in chis do
    try
      Mk := HMFSpace(M, N, k, chi);
      Dk := DihedralBasis(Mk);
    catch e
      Puts(errors_fp, Sprintf("DihedralBasis failed for N=%o (norm %o) chi=%o: %o",
        IdealOneLine(N), Norm(N), HeckeCharLabel(chi : full_label:=false), Sprint(e)));
      Flush(errors_fp);
      continue;
    end try;

    num_dihedral +:= #Dk;

    for f_idx -> f in Dk do
      for bb_idx -> bb in NarrowClassGroupReps(M) do
        coeffs := RestrictionToDiagonal(f, M, bb : AsCoefficients := true);
        stored_coeffs := coeffs[1 .. Min(#coeffs, NUM_COEFFS)];
        is_zero := (#stored_coeffs eq 0) or &and[IsZero(c) : c in stored_coeffs];
        if is_zero then
          num_zero +:= 1;
        end if;

        fname := Sprintf("%o/N%o_chi%o_bb%o_f%o.txt", QEXP_DIR, level_label, chi_idx, bb_idx, f_idx);
        fp := Open(fname, "w");
        Puts(fp, Sprintf("# F = %o, k = %o", F, k));
        Puts(fp, Sprintf("# N = %o (norm %o, label %o)", IdealOneLine(N), Norm(N), level_label));
        Puts(fp, Sprintf("# chi = %o (label %o, order %o)", chi, HeckeCharLabel(chi : full_label:=false), Order(chi)));
        Puts(fp, Sprintf("# bb = %o", IdealOneLine(bb)));
        Puts(fp, Sprintf("# eigenform index (within DihedralBasis) = %o", f_idx));
        Puts(fp, Sprintf("# classical weight = %o, num_coeffs_stored = %o, is_zero = %o",
          n*k[1], #stored_coeffs, is_zero));
        Puts(fp, Join([Sprint(c) : c in stored_coeffs], ","));
        Flush(fp);
      end for;
    end for;
  end for;

  Puts(summary_fp, Sprintf("%o,%o,%o,%o,%o,%o", Norm(N), level_label, IdealOneLine(N), #chis, num_dihedral, num_zero));
  Flush(summary_fp);
  printf "Level %o (norm %o): %o characters, %o dihedral eigenforms, %o zero restrictions\n",
    IdealOneLine(N), Norm(N), #chis, num_dihedral, num_zero;
end for;

print "Done.";
quit;
