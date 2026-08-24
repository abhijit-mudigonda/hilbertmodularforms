load "config.m";

/***************************************************************
* v2: uses the raw-GrpHeckeElt fast-path pipeline (RawDihedralCandidates in
* ModFrmHilD/Creation/Dihedral.m) instead of DihedralBasis/PossibleGrossenchars/
* HMFGrossenchar. This bypasses HMFGrossenchar's Evaluate(chi,I) entirely,
* which is known to crash (silently, if wrapped in try/catch as v1's NewByK
* does) for some (K,D) pairs due to a ClassGroupReps(X)/X`ClassMap desync bug
* -- meaning v1 has been silently under-reporting candidates for any (field,
* level) combination that happens to trigger it. v2 also removes the entire
* outer chi-loop: RawDihedralCandidates finds all psi whose induced
* nebentypus is trivial on (Z/NtgtZ)^* directly, once per (field, target),
* rather than looping over every compatible chi and re-searching per chi.
*
* NEW: curve-coverage tracking. Once every rational newform at a given
* classical level Ntgt has been matched (by ANY field, in this run or a
* previous one), that level is skipped entirely for all subsequent fields --
* there's nothing left to find there. Coverage is tracked in a shared,
* append-only registry file (WeightOneDihedralRestriction/curve_coverage.csv)
* that every parallel worker (sharded by field range) reads fresh before each
* (field, target) pair and appends to on a genuine (scalar != 0) hit. This
* also serves as a persistent index answering "which curves have been fully
* realized" for future runs/queries without needing to re-search.
*
* Output: WeightOneDihedralRestriction/ec_search_v2_m<MIN_M>-<MAX_M>_N<TARGET_LO>-<TARGET_HI>.csv
***************************************************************/

if not assigned TARGET_LO then TARGET_LO := "11"; end if;
if not assigned TARGET_HI then TARGET_HI := "50"; end if;
if not assigned MIN_M then MIN_M := "2"; end if;
if not assigned MAX_M then MAX_M := "100"; end if;
if not assigned DISC_BOUND then DISC_BOUND := "100"; end if;
if not assigned HMF_PRECISION then HMF_PRECISION := "600"; end if;
if not assigned SKIP_PAIRS then SKIP_PAIRS := ""; end if;

TARGET_LO := StringToInteger(TARGET_LO);
TARGET_HI := StringToInteger(TARGET_HI);
MIN_M := StringToInteger(MIN_M);
MAX_M := StringToInteger(MAX_M);
DISC_BOUND := StringToInteger(DISC_BOUND);
HMF_PRECISION := StringToInteger(HMF_PRECISION);

skip_pairs := {};
if SKIP_PAIRS ne "" then
  for tok in Split(SKIP_PAIRS, ",") do
    parts := Split(tok, ":");
    Include(~skip_pairs, <StringToInteger(parts[1]), StringToInteger(parts[2])>);
  end for;
end if;

SetColumns(0);
k11 := [1,1];

OUT_DIR := "WeightOneDihedralRestriction";
System("mkdir -p " * OUT_DIR);
OUT_CSV := Sprintf("%o/ec_search_v2_m%o-%o_N%o-%o.csv", OUT_DIR, MIN_M, MAX_M, TARGET_LO, TARGET_HI);
LOG_FILE := Sprintf("%o/ec_search_v2_m%o-%o_N%o-%o.log", OUT_DIR, MIN_M, MAX_M, TARGET_LO, TARGET_HI);
CHECKPOINT_FILE := Sprintf("%o/ec_search_v2_m%o-%o_N%o-%o.checkpoint", OUT_DIR, MIN_M, MAX_M, TARGET_LO, TARGET_HI);
COVERAGE_FILE := Sprintf("%o/curve_coverage.csv", OUT_DIR);

procedure Log(msg)
  fp := Open(LOG_FILE, "a");
  Puts(fp, msg);
  delete fp;
  print msg;
end procedure;

if OpenTest(CHECKPOINT_FILE, "r") then
  fp := Open(CHECKPOINT_FILE, "r");
  while true do
    line := Gets(fp);
    if IsEof(line) then break; end if;
    parts := Split(line, ",");
    if #parts eq 2 then
      Include(~skip_pairs, <StringToInteger(parts[1]), StringToInteger(parts[2])>);
    end if;
  end while;
  delete fp;
  Log(Sprintf("resumed: %o (m,Ntgt) pairs already done, skipping", #skip_pairs));
end if;

procedure Checkpoint(m, Ntgt)
  fp := Open(CHECKPOINT_FILE, "a");
  Puts(fp, Sprintf("%o,%o", m, Ntgt));
  delete fp;
end procedure;

// -------- curve coverage registry --------
// append-only CSV: target_N,curve_idx,num_curves -- one line per confirmed
// (scalar != 0) hit. Read fresh before each (field,Ntgt) pair so hits found
// by OTHER parallel workers (different field shards) are picked up.
System("test -f '" cat COVERAGE_FILE cat "' || printf 'target_N,curve_idx,num_curves\\n' > '" cat COVERAGE_FILE cat "'");

procedure RecordCoverage(Ntgt, curve_idx, num_curves)
  fp := Open(COVERAGE_FILE, "a");
  Puts(fp, Sprintf("%o,%o,%o", Ntgt, curve_idx, num_curves));
  delete fp;
end procedure;

function IsLevelFullyCovered(Ntgt, num_curves)
  if not OpenTest(COVERAGE_FILE, "r") then
    return false;
  end if;
  matched := {};
  fp := Open(COVERAGE_FILE, "r");
  header_skipped := false;
  while true do
    line := Gets(fp);
    if IsEof(line) then break; end if;
    if not header_skipped then
      header_skipped := true;
      continue;
    end if;
    parts := Split(line, ",");
    if #parts ne 3 then continue; end if;
    if StringToInteger(parts[1]) eq Ntgt then
      Include(~matched, StringToInteger(parts[2]));
    end if;
  end while;
  delete fp;
  return #matched ge num_curves;
end function;

COLS := ["field_disc","target_N","frakN_norm","cand_idx","chi_order",
  "D_norm","mm_norm","num_coeffs","status","ec_dim_at_level","num_rational_newforms",
  "matched_curve_idx","scalar_c","hecke_traces","other_curve_traces",
  "K_rel_disc_norm","K_totally_complex",
  "elapsed_seconds","error_msg"];
NCOLS := #COLS;
C_DISC:=1; C_TARGETN:=2; C_FRAKN_NORM:=3; C_CAND_IDX:=4; C_CHI_ORDER:=5;
C_DNORM:=6; C_MMNORM:=7; C_NUM_COEFFS:=8; C_STATUS:=9; C_ECDIM:=10; C_NUM_RAT:=11;
C_MATCHED_IDX:=12; C_SCALAR_C:=13; C_HECKE_TRACES:=14; C_OTHER_TRACES:=15;
C_K_RELDISC:=16; C_K_TC:=17;
C_ELAPSED:=18; C_ERROR:=19;
if NCOLS ne 19 then error "column mismatch"; end if;

System("test -f '" cat OUT_CSV cat "' || printf '%s\\n' '" cat Join(COLS,",") cat "' > '" cat OUT_CSV cat "'");

function BlankRow()
  return [ "" : i in [1..NCOLS] ];
end function;

procedure WriteRow(row)
  assert #row eq NCOLS;
  fp := Open(OUT_CSV, "a");
  Puts(fp, Join(row, ","));
  delete fp;
end procedure;

function S(x)
  s := Sprint(x);
  s := Join(Split(s, ","), ";");
  s := Join(Split(s, "\n"), " ");
  return s;
end function;

function SerList(seq)
  return Join([Sprint(x) : x in seq], ";");
end function;

ec_cache := AssociativeArray();
procedure GetECSpace(N, ~ec_cache, ~ratnf_out, ~dim_out)
  if IsDefined(ec_cache, N) then
    ratnf_out := ec_cache[N][1];
    dim_out := ec_cache[N][2];
    return;
  end if;
  Mcl := CuspForms(N);
  dim_out := Dimension(Mcl);
  ratnf_out := [* *];
  if dim_out gt 0 then
    nf := Newforms(Mcl);
    for orb in nf do
      g := orb[1];
      if Parent(Coefficient(g,1)) cmpeq RationalField() then
        Append(~ratnf_out, g);
      end if;
    end for;
  end if;
  ec_cache[N] := <ratnf_out, dim_out>;
end procedure;

// ---------------- enumerate real quadratic fields with |disc| <= DISC_BOUND ----------------
field_discs := [];
for m in [MIN_M..MAX_M] do
  if IsSquarefree(m) then
    Fm := QuadraticField(m);
    d := Discriminant(Integers(Fm));
    if d le DISC_BOUND then
      Append(~field_discs, m);
    end if;
  end if;
end for;
Log(Sprintf("=== EC search v2 (raw-GrpHeckeElt fast path, curve-coverage tracking): target levels [%o,%o], %o fields (disc<=%o), precision=%o ===",
  TARGET_LO, TARGET_HI, #field_discs, DISC_BOUND, HMF_PRECISION));
Log(Sprintf("field list (squarefree m): %o", field_discs));

t_start := Cputime();
num_hits := 0;
num_levels_skipped_covered := 0;

targets := [];
for Ntgt in [TARGET_LO..TARGET_HI] do
  GetECSpace(Ntgt, ~ec_cache, ~ratnf0, ~ecdim0);
  if #ratnf0 gt 0 then
    Append(~targets, Ntgt);
  end if;
end for;
Log(Sprintf("of %o candidate levels, %o have >=1 rational newform: %o", TARGET_HI-TARGET_LO+1, #targets, targets));

for m in field_discs do
  t_field := Cputime();
  F := QuadraticField(m);
  ZF := Integers(F);
  fdisc := Discriminant(ZF);
  M := GradedRingOfHMFs(F, HMF_PRECISION);
  // Restricted to the trivial (principal, norm 1) narrow class component only.
  // RestrictionToDiagonal's "exp := j div b" for a nontrivial component bb
  // (b := generator of bb meet Z) incorrectly collapses coefficients of what
  // is actually a level-(b*Ntgt) form into a mislabeled level-Ntgt one --
  // b=1 exactly when bb=ZF, so only the trivial component is currently
  // trustworthy. Skipping nontrivial components entirely for now rather
  // than risk false or corrupted matches; revisit once RestrictionToDiagonal
  // is fixed to handle bb!=1 correctly.
  bbs := [bb : bb in NarrowClassGroupReps(M) | Norm(bb) eq 1];
  assert #bbs eq 1;

  for Ntgt in targets do
    if <m, Ntgt> in skip_pairs then
      Log(Sprintf("    field Q(sqrt%o) target_N=%o: SKIPPED (in SKIP_PAIRS/checkpoint)", m, Ntgt));
      continue;
    end if;
    t_tgt := Cputime();
    GetECSpace(Ntgt, ~ec_cache, ~ratnf, ~ecdim);

    if IsLevelFullyCovered(Ntgt, #ratnf) then
      num_levels_skipped_covered +:= 1;
      Log(Sprintf("    field Q(sqrt%o) target_N=%o: SKIPPED (all %o curve(s) already covered)", m, Ntgt, #ratnf));
      Checkpoint(m, Ntgt);
      continue;
    end if;

    frakN := ideal<ZF | Ntgt>;

    try
      Dk := RawDihedralCandidates(M, frakN, Ntgt);
    catch e
      Log(Sprintf("    field Q(sqrt%o) target_N=%o: ERROR in RawDihedralCandidates: %o", m, Ntgt, e));
      Checkpoint(m, Ntgt);
      continue;
    end try;

    if #Dk eq 0 then
      Checkpoint(m, Ntgt);
      continue;
    end if;

    base := BlankRow();
    base[C_DISC] := S(fdisc);
    base[C_TARGETN] := S(Ntgt);
    base[C_FRAKN_NORM] := S(Norm(frakN));
    base[C_ECDIM] := S(ecdim);
    base[C_NUM_RAT] := S(#ratnf);

    for cand_idx -> cand in Dk do
      f, D, K, mm := Explode(cand);
      row := base;
      row[C_CAND_IDX] := S(cand_idx);
      row[C_DNORM] := S(Norm(D));
      row[C_MMNORM] := S(Norm(mm));
      row[C_K_RELDISC] := S(Norm(Discriminant(Integers(K))));
      row[C_K_TC] := S(IsTotallyComplex(K));
      try
        row[C_CHI_ORDER] := S(Order(Character(Parent(f))));
      catch e
        row[C_CHI_ORDER] := "?";
      end try;

      for bb_idx -> bb in bbs do
        t1 := Cputime();
        rowb := row;

        try
          coeffs := RestrictionToDiagonal(f, M, bb : AsCoefficients := true);
        catch e
          rowb[C_STATUS] := "err_restriction";
          rowb[C_ELAPSED] := S(Cputime()-t1);
          rowb[C_ERROR] := S(e);
          WriteRow(rowb);
          continue;
        end try;

        if forall{c : c in coeffs | c eq 0} then
          rowb[C_STATUS] := "zero_restriction";
          rowb[C_ELAPSED] := S(Cputime()-t1);
          WriteRow(rowb);
          continue;
        end if;

        is_rat := true;
        for c in coeffs do
          if not IsCoercible(RationalField(), c) then is_rat := false; break; end if;
        end for;
        if not is_rat then
          rowb[C_STATUS] := "non_rational_coeffs";
          rowb[C_ELAPSED] := S(Cputime()-t1);
          WriteRow(rowb);
          continue;
        end if;

        target := [RationalField()!c : c in coeffs];
        NMAX := #target - 1;
        v := [target[idx+1] : idx in [1..NMAX]];
        rowb[C_NUM_COEFFS] := S(#target);

        // A degeneracy lift (mm != 1) forces a_n(restriction) = 0 for every
        // n not divisible by (essentially) Norm(mm) -- a_1..a_{Norm(mm)-1}
        // are structurally zero regardless of the underlying form, so they
        // carry no information. Requiring NMAX > ecdim alone is not enough:
        // it can pass while v is still entirely within that structurally-zero
        // range, producing a spurious scalar=0 "match" via the trivial
        // solution to IsConsistent. Require real precision headroom past the
        // degeneracy shift, with a margin so the match rests on several
        // genuinely-informative (potentially nonzero) coefficients, not just
        // clearing a bare threshold.
        min_needed := Max(ecdim, Norm(mm)) + 5;
        if NMAX le min_needed then
          rowb[C_STATUS] := "inconclusive_insufficient_precision";
          rowb[C_ELAPSED] := S(Cputime()-t1);
          WriteRow(rowb);
          continue;
        end if;

        if forall{c : c in v | c eq 0} then
          rowb[C_STATUS] := "inconclusive_zero_after_precision";
          rowb[C_ELAPSED] := S(Cputime()-t1);
          WriteRow(rowb);
          continue;
        end if;

        matched := false;
        for ci -> g in ratnf do
          a := [Rationals()!Coefficient(g, idx) : idx in [1..NMAX]];
          T := Matrix(Rationals(), 1, NMAX, a);
          ok, csol := IsConsistent(T, Vector(Rationals(), v));
          // scalar=0 is a degenerate/spurious "match" (only possible when v
          // itself is the zero vector, i.e. every nonzero Hecke eigenvalue
          // in the restriction is zero) -- not a genuine realization, so it
          // does not count and must not be recorded as coverage. (Guarded
          // above too, but kept here as a belt-and-suspenders check.)
          if ok and csol[1] ne 0 then
            c := csol[1];
            rowb[C_STATUS] := "matched_elliptic_curve";
            rowb[C_MATCHED_IDX] := S(ci);
            rowb[C_SCALAR_C] := S(c);
            small_primes := [p : p in PrimesUpTo(40) | Ntgt mod p ne 0];
            rowb[C_HECKE_TRACES] := SerList([Coefficient(g,p) : p in small_primes]);
            if #ratnf gt 1 then
              other := [<cj, SerList([Coefficient(ratnf[cj],p) : p in small_primes])> : cj in [1..#ratnf] | cj ne ci];
              rowb[C_OTHER_TRACES] := S(other);
            end if;
            matched := true;
            num_hits +:= 1;
            RecordCoverage(Ntgt, ci, #ratnf);

            Log(Sprintf(">>> HIT: field disc=%o target_N=%o matched curve #%o of %o, scalar=%o, D_norm=%o mm_norm=%o K_reldisc=%o totally_complex=%o",
              fdisc, Ntgt, ci, #ratnf, c, Norm(D), Norm(mm), Norm(Discriminant(Integers(K))), IsTotallyComplex(K)));
            break;
          end if;
        end for;
        if not matched then
          rowb[C_STATUS] := "not_ec_shape";
        end if;
        rowb[C_ELAPSED] := S(Cputime()-t1);
        WriteRow(rowb);
      end for;
    end for;
    Log(Sprintf("    field Q(sqrt%o) target_N=%o: %o candidates, elapsed %o s (field total %o s)",
      m, Ntgt, #Dk, Cputime()-t_tgt, Cputime()-t_field));
    Checkpoint(m, Ntgt);
  end for;
  Log(Sprintf("done field Q(sqrt%o) [disc %o], elapsed %o s (total %o s), hits so far %o, levels skipped as fully-covered %o",
    m, fdisc, Cputime()-t_field, Cputime()-t_start, num_hits, num_levels_skipped_covered));
end for;

Log(Sprintf("=== EC search v2 complete: %o hits, %o levels skipped as fully-covered, total elapsed %o s ===",
  num_hits, num_levels_skipped_covered, Cputime()-t_start));
quit;
