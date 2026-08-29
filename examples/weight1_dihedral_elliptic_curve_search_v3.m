load "config.m";

/***************************************************************
* v3: uses the ideal-coefficient bypass pipeline (RawDihedralCurveMatch in
* ModFrmHilD/Creation/Dihedral.m) instead of GradedRingOfHMFs +
* RawDihedralCandidates + RestrictionToDiagonal + IsConsistent. The bypass
* computes only the trace-indexed classical coefficients the diagonal
* restriction actually needs (no ambient HMF norm-precision, no
* GradedRingOfHMFs at all), with early exit as soon as a candidate provably
* can't match any target curve. This also removes the crash risk v2 had for
* fields with large fundamental units, where HMFConstructExpansion's
* Monomial(R,e) storage could exceed Magma's exponent ceiling even for
* mathematically ordinary elements.
*
* JMAX (how many classical Fourier coefficients to require agreement on
* before trusting a match) is chosen per target level automatically, as
* Dimension(CuspForms(Ntgt)) + a margin -- mirroring v2's own
* min_needed := Max(ecdim, Norm(mm)) + 5 precision-sufficiency threshold.
* Phase 1+2+3 are all in scope: direct D=frakN matches and Inclusion-style
* oldform lifts for every divisor dd | frakN/D, narrowly principal or not
* (empirically validated against the real Inclusion intrinsic and the full
* historical hit set -- see session notes).
*
* Curve-coverage tracking (same mechanism as v2): once every rational
* newform at a given classical level Ntgt has been matched (by ANY field, in
* this run or a previous one), that level is skipped entirely for all
* subsequent fields. Coverage is tracked in a shared, append-only registry
* file (WeightOneDihedralRestriction/curve_coverage.csv) that every parallel
* worker (sharded by field range) reads fresh before each (field, target)
* pair and appends to on a genuine (scalar != 0) hit.
*
* Output: WeightOneDihedralRestriction/ec_search_v3_m<MIN_M>-<MAX_M>_N<TARGET_LO>-<TARGET_HI>.csv
***************************************************************/

if not assigned TARGET_LO then TARGET_LO := "11"; end if;
if not assigned TARGET_HI then TARGET_HI := "50"; end if;
if not assigned MIN_M then MIN_M := "2"; end if;
if not assigned MAX_M then MAX_M := "100"; end if;
if not assigned DISC_BOUND then DISC_BOUND := "100"; end if;
if not assigned JMAX_MARGIN then JMAX_MARGIN := "10"; end if;
if not assigned SKIP_PAIRS then SKIP_PAIRS := ""; end if;
if not assigned M_STRIDE then M_STRIDE := "1"; end if;
if not assigned M_RESIDUE then M_RESIDUE := "0"; end if;

TARGET_LO := StringToInteger(TARGET_LO);
TARGET_HI := StringToInteger(TARGET_HI);
MIN_M := StringToInteger(MIN_M);
MAX_M := StringToInteger(MAX_M);
DISC_BOUND := StringToInteger(DISC_BOUND);
JMAX_MARGIN := StringToInteger(JMAX_MARGIN);
M_STRIDE := StringToInteger(M_STRIDE);
M_RESIDUE := StringToInteger(M_RESIDUE);

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
OUT_TAG := (M_STRIDE gt 1) select Sprintf("m%o-%o_r%oof%o", MIN_M, MAX_M, M_RESIDUE, M_STRIDE) else Sprintf("m%o-%o", MIN_M, MAX_M);
OUT_CSV := Sprintf("%o/ec_search_v3_%o_N%o-%o.csv", OUT_DIR, OUT_TAG, TARGET_LO, TARGET_HI);
LOG_FILE := Sprintf("%o/ec_search_v3_%o_N%o-%o.log", OUT_DIR, OUT_TAG, TARGET_LO, TARGET_HI);
CHECKPOINT_FILE := Sprintf("%o/ec_search_v3_%o_N%o-%o.checkpoint", OUT_DIR, OUT_TAG, TARGET_LO, TARGET_HI);
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

// -------- curve coverage registry (shared with v2) --------
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

COLS := ["field_disc","target_N","frakN_norm","ec_dim_at_level","num_rational_newforms",
  "JMAX","K_rel_disc_norm","K_totally_complex","D_norm","dd_norm",
  "matched_curve_idx","scalar_c","died_at","hecke_traces","other_curve_traces",
  "elapsed_seconds"];
NCOLS := #COLS;
C_DISC:=1; C_TARGETN:=2; C_FRAKN_NORM:=3; C_ECDIM:=4; C_NUM_RAT:=5;
C_JMAX:=6; C_K_RELDISC:=7; C_K_TC:=8; C_DNORM:=9; C_MMNORM:=10;
C_MATCHED_IDX:=11; C_SCALAR_C:=12; C_DIED_AT:=13; C_HECKE_TRACES:=14; C_OTHER_TRACES:=15;
C_ELAPSED:=16;
if NCOLS ne 16 then error "column mismatch"; end if;

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
all_field_discs := [];
for m in [MIN_M..MAX_M] do
  if IsSquarefree(m) then
    Fm := QuadraticField(m);
    d := Discriminant(Integers(Fm));
    if d le DISC_BOUND then
      Append(~all_field_discs, m);
    end if;
  end if;
end for;
field_discs := [all_field_discs[i] : i in [1..#all_field_discs] | ((i-1) mod M_STRIDE) eq M_RESIDUE];
Log(Sprintf("=== EC search v3 (ideal-coefficient bypass, curve-coverage tracking): target levels [%o,%o], %o fields (disc<=%o), JMAX_MARGIN=%o ===",
  TARGET_LO, TARGET_HI, #field_discs, DISC_BOUND, JMAX_MARGIN));
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
    JMAX := ecdim + JMAX_MARGIN;

    try
      matches := RawDihedralCurveMatch(F, frakN, Ntgt, ratnf : JMAX := JMAX);
    catch e
      Log(Sprintf("    field Q(sqrt%o) target_N=%o: ERROR in RawDihedralCurveMatch: %o", m, Ntgt, e));
      Checkpoint(m, Ntgt);
      continue;
    end try;

    for mtch in matches do
      K, Dnorm, ddnorm, ci, scalar, died_at, outcome := Explode(mtch);
      if scalar eq 0 then continue; end if; // spurious/degenerate, should not occur but guard anyway

      row := BlankRow();
      row[C_DISC] := S(fdisc);
      row[C_TARGETN] := S(Ntgt);
      row[C_FRAKN_NORM] := S(Norm(frakN));
      row[C_ECDIM] := S(ecdim);
      row[C_NUM_RAT] := S(#ratnf);
      row[C_JMAX] := S(JMAX);
      row[C_K_RELDISC] := S(Norm(Discriminant(Integers(K))));
      row[C_K_TC] := S(IsTotallyComplex(K));
      row[C_DNORM] := S(Dnorm);
      row[C_MMNORM] := S(ddnorm);
      row[C_MATCHED_IDX] := S(ci);
      row[C_SCALAR_C] := S(scalar);
      row[C_DIED_AT] := S(died_at);

      small_primes := [p : p in PrimesUpTo(40) | Ntgt mod p ne 0];
      row[C_HECKE_TRACES] := SerList([Coefficient(ratnf[ci],p) : p in small_primes]);
      if #ratnf gt 1 then
        other := [<cj, SerList([Coefficient(ratnf[cj],p) : p in small_primes])> : cj in [1..#ratnf] | cj ne ci];
        row[C_OTHER_TRACES] := S(other);
      end if;
      row[C_ELAPSED] := S(Cputime()-t_tgt);
      WriteRow(row);

      num_hits +:= 1;
      RecordCoverage(Ntgt, ci, #ratnf);
      Log(Sprintf(">>> HIT: field disc=%o target_N=%o matched curve #%o of %o, scalar=%o, D_norm=%o dd_norm=%o K_reldisc=%o totally_complex=%o JMAX=%o",
        fdisc, Ntgt, ci, #ratnf, scalar, Dnorm, ddnorm, Norm(Discriminant(Integers(K))), IsTotallyComplex(K), JMAX));
    end for;

    Log(Sprintf("    field Q(sqrt%o) target_N=%o: %o matches, JMAX=%o, elapsed %o s (field total %o s)",
      m, Ntgt, #matches, JMAX, Cputime()-t_tgt, Cputime()-t_field));
    Checkpoint(m, Ntgt);
  end for;
  Log(Sprintf("done field Q(sqrt%o) [disc %o], elapsed %o s (total %o s), hits so far %o, levels skipped as fully-covered %o",
    m, fdisc, Cputime()-t_field, Cputime()-t_start, num_hits, num_levels_skipped_covered));
end for;

Log(Sprintf("=== EC search v3 complete: %o hits, %o levels skipped as fully-covered, total elapsed %o s ===",
  num_hits, num_levels_skipped_covered, Cputime()-t_start));
quit;
