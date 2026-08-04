load "config.m";

/***************************************************************
* Surveys parallel weight [1,1] dihedral Hilbert modular eigenforms
* over a real quadratic field F, across levels N up to a norm bound,
* restricts each eigenform to the diagonal, and attempts to identify
* the restriction as Tr_{K/Q}(c*g) for a classical weight 2 newform
* orbit g (K = coefficient field of g).
*
* Usage (parameters are passed on the command line and read back as
* strings, so we StringToInteger them below):
*
*   magma FIELD_DISC:=5 LEVEL_BOUND:=1500 examples/weight1_dihedral_classical_survey.m
*
*   FIELD_DISC     -- 5 for Q(sqrt5), 2 for Q(sqrt2)
*   LEVEL_BOUND    -- upper bound on Norm(N)
*   NORM_MIN       -- lower bound on Norm(N) (default 1); use this to
*                     split a run into ranges across parallel jobs
*   HMF_PRECISION  -- OPTIONAL. If set to an integer, disables adaptive
*                     precision and uses this single fixed Fourier
*                     precision for every level (the old behaviour --
*                     useful for a targeted rerun of one troublesome
*                     range). If unset (default), precision is chosen
*                     ADAPTIVELY per level -- see below.
*   MIN_PRECISION  -- floor on the adaptive precision (default 300).
*   MAX_PRECISION  -- ceiling on the adaptive precision (default 60000);
*                     levels whose dimension-based estimate exceeds this
*                     are processed at the ceiling anyway and are likely
*                     to come back "inconclusive_insufficient_precision"
*                     -- rerun those specific levels by hand with
*                     HMF_PRECISION set explicitly if you care about them.
*
* ADAPTIVE PRECISION: the ambient graded ring's Fourier precision is the
* single biggest cost lever in this script, and the number of restriction
* coefficients obtained grows only like sqrt(precision) (see
* weight1_dihedral_restriction_Qsqrt5.m's comment), so precision needs to
* grow QUADRATICALLY with the number of coefficients we want -- and the
* number of coefficients we need scales with the dimension of the
* classical cusp space we're testing membership in. Using one large fixed
* precision for every level (the previous version of this script) wastes
* most of its time on low-norm levels whose classical space has dimension
* 1-5. Instead, for each ideal N we:
*   1. cheaply compute Dimension(CuspForms(L)) and Dimension of a
*      representative even-character space for each nebentypus order
*      that will actually be searched (these are closed-form dimension
*      formula lookups, NOT Newforms computations, so this is fast even
*      up to level 1500);
*   2. convert the largest such dimension into a target coefficient count
*      (with safety margin) and invert the empirical sqrt scaling law to
*      get a needed precision, rounded up to a shared power-of-two tier
*      (300, 600, 1200, 2400, ...);
*   3. bucket all ideals by their assigned tier and build ONE graded ring
*      per tier (not per ideal -- rebuilding the ring itself has fixed
*      overhead, so we still want to amortize it over many ideals), then
*      process each tier's ideals against that tier's ring.
* This means ideals are processed in tier order (cheap tiers first, so
* low-norm/low-dimension levels finish almost immediately), not strictly
* in norm order -- the per-level log lines will jump around in norm
* within a tier's block, which is expected.
*
* Run the two fields in separate processes to parallelize:
*   magma FIELD_DISC:=5 LEVEL_BOUND:=1500 examples/weight1_dihedral_classical_survey.m &
*   magma FIELD_DISC:=2 LEVEL_BOUND:=1500 examples/weight1_dihedral_classical_survey.m &
*
* Output (written incrementally, one flushed row at a time, so a
* killed/crashed run keeps all progress up to that point):
*   WeightOneDihedralRestriction/survey_Qsqrt<D>_<NORM_MIN>-<LEVEL_BOUND>.csv
*   WeightOneDihedralRestriction/survey_Qsqrt<D>_<NORM_MIN>-<LEVEL_BOUND>.log
*
* IMPORTANT CAVEATS (read before trusting a row):
*  - "candidate_label" is NOT a verified LMFDB label. It is a local
*    identifier of the form "<L>.2.o<chi_order>i<chi_idx>.orb<orbit_idx>"
*    encoding enough data (level, weight 2, nebentypus order/index,
*    orbit index) to help a human re-find the orbit, but the exact
*    LMFDB orbit letter requires either network access or
*    reimplementing LMFDB's canonical ordering, neither of which this
*    script does. Use minpoly_K + hecke_traces to look it up by hand
*    (or have Claude do it) after the fact.
*  - The classical nebentypus search only tries characters whose order
*    divides Order(chi) (the HMF nebentypus order), plus the trivial
*    character. This is a heuristic motivated by every hand-checked
*    example so far (chi order 2 <-> classical nebentypus order 2 or
*    1), NOT a theorem. Rows with status "unresolved_tried_divisors"
*    are cases where this heuristic failed to find a match -- they are
*    not necessarily unidentifiable, just not identified by this
*    search. Cross check by hand if there are many of these.
*  - "matched" only means IsConsistent found an EXACT rational-linear-
*    algebra fit AND num_coeffs strictly exceeds the classical space's
*    dimension (i.e. the test is overdetermined, not tautological).
*    If num_coeffs <= classical_dim for every candidate tried, status is
*    "inconclusive_insufficient_precision" -- the adaptive precision
*    estimate undershot for that level (its true classical dimension
*    exceeded what the cheap dimension-formula-based estimate predicted,
*    e.g. because the actual classical level differs from the Norm(N)
*    proxy used to size it). Rerun just that norm range with
*    HMF_PRECISION set explicitly higher to resolve it.
*  - The (M, psi) columns (the CM/RM quadratic extension of F and the
*    Grossencharacter used to build the theta series) are matched
*    best-effort by comparing Hecke eigenvalues of f against
*    ThetaSeries(Mk, psi) for psi in PossibleGrossenchars(Mk) at a
*    handful of primes. This only searches characters at the CURRENT
*    level N; eigenforms coming from OldDihedralBasis (lifted from a
*    proper divisor of N) will generally not match anything here and
*    are logged as such -- this is expected, not a bug.
***************************************************************/

if not assigned FIELD_DISC then FIELD_DISC := "5"; end if;
if not assigned LEVEL_BOUND then LEVEL_BOUND := "1500"; end if;
if not assigned NORM_MIN then NORM_MIN := "1"; end if;
if not assigned HMF_PRECISION then HMF_PRECISION := "auto"; end if;
if not assigned MIN_PRECISION then MIN_PRECISION := "300"; end if;
if not assigned MAX_PRECISION then MAX_PRECISION := "60000"; end if;

D := StringToInteger(FIELD_DISC);
LEVEL_BOUND := StringToInteger(LEVEL_BOUND);
NORM_MIN := StringToInteger(NORM_MIN);
MIN_PRECISION := StringToInteger(MIN_PRECISION);
MAX_PRECISION := StringToInteger(MAX_PRECISION);

ADAPTIVE_PRECISION := (HMF_PRECISION eq "auto");
// NB: FIXED_PRECISION must be defined unconditionally (sentinel 0 when
// unused) rather than only inside "if not ADAPTIVE_PRECISION". Magma
// appears to validate identifier bindings across ALL branches of an
// if/else loaded from a file -- even the branch that won't run -- and
// silently aborts the whole if/else statement if any branch references
// an identifier that isn't unconditionally defined. Leaving this
// conditionally-defined caused every later "if ADAPTIVE_PRECISION then
// ... else ... (uses FIXED_PRECISION) ... end if;" to fail outright.
FIXED_PRECISION := ADAPTIVE_PRECISION select 0 else StringToInteger(HMF_PRECISION);

if not (D eq 2 or D eq 5) then
  error "FIELD_DISC must be 2 or 5 (this script is scoped to Q(sqrt2), Q(sqrt5))";
end if;

SetColumns(0);

F := QuadraticField(D);
ZF := Integers(F);
k := [1, 1];
n := Degree(F);

field_tag := Sprintf("Qsqrt%o", D);
OUT_DIR := "WeightOneDihedralRestriction";
System("mkdir -p " * OUT_DIR);
OUT_CSV := Sprintf("%o/survey_%o_%o-%o.csv", OUT_DIR, field_tag, NORM_MIN, LEVEL_BOUND);
LOG_FILE := Sprintf("%o/survey_%o_%o-%o.log", OUT_DIR, field_tag, NORM_MIN, LEVEL_BOUND);

// Column layout. WriteRow always receives a full-length blank-initialized
// row (see BlankRow) with individual slots set by column index -- this
// avoids the (very real, previously present) bug class of miscounting
// commas in giant hand-written literal row tuples.
COLS := ["field","hmf_level_label","hmf_level_norm","chi_idx","num_chis","chi_order","chi_label","f_idx","num_f","bb_idx",
  "num_coeffs","coeffs_integral","classical_level","status","neb_order","neb_char_idx","classical_dim","num_orbits",
  "matched_orbit_idx","deg_K","minpoly_K","num_aut_K","minpoly_c","deg_c","trace_c",
  "psi_match_found","psi_order","psi_conductor_norm","M_defining_poly","M_disc","M_totally_complex","M_class_number",
  "candidate_label","hecke_traces","elapsed_seconds","error_msg"];
NCOLS := #COLS;
HEADER := Join(COLS, ",");

C_FIELD:=1; C_LEVEL_LABEL:=2; C_LEVEL_NORM:=3; C_CHI_IDX:=4; C_NUM_CHIS:=5; C_CHI_ORDER:=6; C_CHI_LABEL:=7;
C_F_IDX:=8; C_NUM_F:=9; C_BB_IDX:=10; C_NUM_COEFFS:=11; C_COEFFS_INT:=12; C_CLASSICAL_LEVEL:=13; C_STATUS:=14;
C_NEB_ORDER:=15; C_NEB_CHAR_IDX:=16; C_CLASSICAL_DIM:=17; C_NUM_ORBITS:=18; C_MATCHED_ORBIT_IDX:=19; C_DEG_K:=20;
C_MINPOLY_K:=21; C_NUM_AUT_K:=22; C_MINPOLY_C:=23; C_DEG_C:=24; C_TRACE_C:=25; C_PSI_FOUND:=26; C_PSI_ORDER:=27;
C_PSI_COND_NORM:=28; C_M_POLY:=29; C_M_DISC:=30; C_M_TC:=31; C_M_CLASSNUM:=32; C_CANDIDATE_LABEL:=33;
C_HECKE_TRACES:=34; C_ELAPSED:=35; C_ERROR:=36;

if NCOLS ne 36 then
  error "column count / index constants are out of sync -- fix COLS or the C_* constants";
end if;

System("test -f '" cat OUT_CSV cat "' || printf '%s\\n' '" cat HEADER cat "' > '" cat OUT_CSV cat "'");

procedure Log(msg)
  fp := Open(LOG_FILE, "a");
  Puts(fp, msg);
  delete fp;
  print msg;
end procedure;

function BlankRow()
  return [ "" : i in [1..NCOLS] ];
end function;

procedure WriteRow(row)
  // row: SeqEnum of length NCOLS of already-stringified column values
  assert #row eq NCOLS;
  fp := Open(OUT_CSV, "a");
  Puts(fp, Join(row, ","));
  delete fp;
end procedure;

function S(x)
  // stringify + neutralize commas/newlines so the CSV stays well-formed
  s := Sprint(x);
  s := Join(Split(s, ","), ";");
  s := Join(Split(s, "\n"), " ");
  return s;
end function;

function SerList(seq)
  return Join([Sprint(x) : x in seq], ";");
end function;

// ---------------- caches, keyed by string ----------------
// Magma functions cannot silently mutate an outer-scope variable, so the
// cache is threaded through explicitly as a reference (~) parameter.
classical_cache := AssociativeArray(); // key -> <nf, dim>

procedure GetClassicalSpace(L, order, char_idx_in_evens, ~classical_cache, ~nf_out, ~dim_out)
  // order eq 1 means trivial. Otherwise char_idx_in_evens indexes into
  // the list of even Dirichlet characters mod L of order exactly `order`.
  // Magma only allows ~ (reference) parameters on procedures, not
  // functions, hence this is a procedure with output parameters rather
  // than a function -- it needs write access to classical_cache.
  key := Sprintf("%o_%o_%o", L, order, char_idx_in_evens);
  if IsDefined(classical_cache, key) then
    nf_out := classical_cache[key][1];
    dim_out := classical_cache[key][2];
    return;
  end if;
  if order eq 1 then
    Mcl := CuspForms(L);
  else
    G := DirichletGroup(L);
    evens := [c : c in Elements(G) | Order(c) eq order and Evaluate(c,-1) eq 1];
    Mcl := CuspidalSubspace(ModularForms(evens[char_idx_in_evens], 2));
  end if;
  nf_out := Newforms(Mcl);
  dim_out := Dimension(Mcl);
  classical_cache[key] := <nf_out, dim_out>;
end procedure;

function NumEvenCharsOfOrder(L, order)
  if order eq 1 then return 1; end if;
  G := DirichletGroup(L);
  return #[c : c in Elements(G) | Order(c) eq order and Evaluate(c,-1) eq 1];
end function;

// ---------------- identification core ----------------

function AnalyzeMembership(nf, v, NMAX)
  // Returns: ok, dim, num_orbits, contributing_orbit_idx (0 if multiple/none),
  // deg_K, minpoly_K_string, num_aut_K, c_elt, minpoly_c_string, deg_c, trace_c
  cols := [];
  for i in [1..#nf] do
    f := nf[i][1];
    K := Parent(Coefficient(f,1));
    d := Degree(K);
    a := [Coefficient(f,idx) : idx in [1..NMAX]];
    if d eq 1 then
      Append(~cols, [Rationals()!a[idx] : idx in [1..NMAX]]);
    else
      for j in [0..d-1] do
        Append(~cols, [Rationals()!Trace(K.1^j * a[idx]) : idx in [1..NMAX]]);
      end for;
    end if;
  end for;
  dim := #cols;
  if dim eq 0 or NMAX le dim then
    return false, dim, #nf, 0, 0, "", 0, "", "", 0, "";
  end if;
  Tfull := Matrix(Rationals(), cols);
  ok, csol := IsConsistent(Tfull, Vector(Rationals(), v));
  if not ok then
    return false, dim, #nf, 0, 0, "", 0, "", "", 0, "";
  end if;
  // which orbit(s) contribute
  pos := 1;
  contributing := [];
  orbit_c := AssociativeArray();
  for i in [1..#nf] do
    K := Parent(Coefficient(nf[i][1],1));
    d := Degree(K);
    block := [csol[pos+j] : j in [0..d-1]];
    if block ne [Rationals()|0 : x in block] then
      Append(~contributing, i);
      orbit_c[i] := <K, block, d>;
    end if;
    pos +:= d;
  end for;
  if #contributing ne 1 then
    // matches only as a combination across multiple orbits; still log it
    // but without a single clean (K, c)
    return true, dim, #nf, 0, 0, "", 0, "", "", 0, SerList(contributing);
  end if;
  orbit_idx := contributing[1];
  K, block, d := Explode(orbit_c[orbit_idx]);
  if d eq 1 then
    c := block[1];
    deg_c := 1;
    minpoly_c_str := Sprint(MinimalPolynomial(c));
    trace_c := c;
    naut := 1;
  else
    c := &+[block[j+1]*K.1^j : j in [0..d-1]];
    mpc := MinimalPolynomial(c);
    deg_c := Degree(mpc);
    minpoly_c_str := Sprint(mpc);
    trace_c := Trace(c);
    naut := #Automorphisms(K);
  end if;
  minpoly_K_str := (d eq 1) select "x" else Sprint(MinimalPolynomial(K.1));
  return true, dim, #nf, orbit_idx, d, minpoly_K_str, naut, minpoly_c_str, deg_c, trace_c, "";
end function;

function DivisorsOf(m)
  return Sort(Divisors(m));
end function;

procedure TryIdentify(L, chi_order, v, NMAX, ~classical_cache,
    ~status, ~neb_order, ~neb_char_idx, ~classical_dim, ~num_orbits,
    ~matched_orbit_idx, ~minpoly_K, ~deg_K, ~minpoly_c, ~contrib, ~deg_c, ~trace_c, ~num_aut_K)
  divs := DivisorsOf(chi_order);
  CANDIDATE_CAP := 40; // soft cap on number of distinct classical characters tried
  tried := 0;
  best_inconclusive := false;
  for ord in divs do
    nchars := NumEvenCharsOfOrder(L, ord);
    for cidx in [1..nchars] do
      tried +:= 1;
      if tried gt CANDIDATE_CAP then
        status:="unresolved_hit_candidate_cap"; neb_order:=ord; neb_char_idx:=cidx;
        classical_dim:=0; num_orbits:=0; matched_orbit_idx:=0; minpoly_K:="";
        deg_K:=0; minpoly_c:=""; contrib:=""; deg_c:=0; trace_c:=0; num_aut_K:=0;
        return;
      end if;
      GetClassicalSpace(L, ord, cidx, ~classical_cache, ~nf, ~dim);
      ok, dimr, norb, orbit_idx, degK, minpolyK, naut, minpolyc, degc, tracec, contrib_local :=
        AnalyzeMembership(nf, v, NMAX);
      if NMAX le dimr then
        best_inconclusive := true;
      end if;
      if ok then
        if orbit_idx gt 0 then
          status:="matched"; neb_order:=ord; neb_char_idx:=cidx; classical_dim:=dimr;
          num_orbits:=norb; matched_orbit_idx:=orbit_idx; minpoly_K:=minpolyK; deg_K:=degK;
          minpoly_c:=minpolyc; contrib:=""; deg_c:=degc; trace_c:=tracec; num_aut_K:=naut;
          return;
        else
          status:="matched_multi_orbit"; neb_order:=ord; neb_char_idx:=cidx; classical_dim:=dimr;
          num_orbits:=norb; matched_orbit_idx:=0; minpoly_K:=""; deg_K:=0; minpoly_c:="";
          contrib:=contrib_local; deg_c:=0; trace_c:=0; num_aut_K:=0;
          return;
        end if;
      end if;
    end for;
  end for;
  status := best_inconclusive select "inconclusive_insufficient_precision" else "unresolved_tried_divisors";
  neb_order:=0; neb_char_idx:=0; classical_dim:=0; num_orbits:=0; matched_orbit_idx:=0;
  minpoly_K:=""; deg_K:=0; minpoly_c:=""; contrib:=""; deg_c:=0; trace_c:=0; num_aut_K:=0;
end procedure;

// ---------------- (M, psi) best-effort matching ----------------

function TryMatchGrossenchar(Mk, f)
  try
    psis := PossibleGrossenchars(Mk);
  catch e
    return false, 0, 0, "", 0, false, 0;
  end try;
  if #psis eq 0 then
    return false, 0, 0, "", 0, false, 0;
  end if;
  pps := [pp : pp in PrimeIdeals(Parent(Mk)) | Norm(pp) le 80 and Norm(pp) gt 1];
  if #pps eq 0 then
    return false, 0, 0, "", 0, false, 0;
  end if;
  fcoeffs := [Coefficient(f, pp) : pp in pps];
  for psi in psis do
    try
      g := ThetaSeries(Mk, psi);
      gcoeffs := [Coefficient(g, pp) : pp in pps];
    catch e
      continue;
    end try;
    if gcoeffs eq fcoeffs then
      Kcm := AbsoluteField(NumberField(Order(Modulus(psi))));
      psi_ord := Order(psi`RayClassChar);
      cond := Conductor(psi);
      cond_norm := Norm(cond);
      is_tc := IsTotallyComplex(Kcm);
      class_num := 0;
      if cond_norm eq 1 then
        try
          ClM := ClassGroup(Kcm);
          class_num := #ClM;
        catch e
          class_num := -1; // computation failed/too slow, sentinel
        end try;
      end if;
      return true, psi_ord, cond_norm, Sprint(DefiningPolynomial(Kcm)), Discriminant(Integers(Kcm)), is_tc, class_num;
    end if;
  end for;
  return false, 0, 0, "", 0, false, 0;
end function;

// ---------------- adaptive precision planning ----------------
// See the file header comment for the rationale. All of this only uses
// closed-form Dimension() lookups (no Newforms/eigenform computation),
// so it's cheap even run once per ideal across the whole norm range.

// empirical: coeffs(precision) ~ 0.894*sqrt(precision), i.e.
// precision ~ (coeffs/0.894)^2 ~ 1.251*coeffs^2. We inflate the constant
// a bit (1.4 instead of 1.251) since this is a rough fit from a handful
// of data points, not a proven law.
PRECISION_SCALE_NUM := 7; PRECISION_SCALE_DEN := 5; // 1.4 as a rational
MARGIN_MULT_NUM := 13; MARGIN_MULT_DEN := 10; // 1.3x the raw dimension estimate
BASE_MARGIN_ADD := 10; // plus an additive safety buffer

function EstimateClassicalDim(L, chi_orders)
  // chi_orders: SetEnum/SeqEnum of the HMF nebentypus orders that will
  // actually be searched for this ideal (see TryIdentify's divisor-of-
  // chi-order heuristic) -- we size precision off the largest classical
  // space among {trivial} union {one representative even character of
  // each of these orders, and of each of their divisors}.
  best := Dimension(CuspForms(L));
  orders_to_check := {1};
  for co in chi_orders do
    for d in DivisorsOf(co) do
      Include(~orders_to_check, d);
    end for;
  end for;
  G := DirichletGroup(L);
  evens_by_order := AssociativeArray();
  for chi in Elements(G) do
    if Evaluate(chi,-1) eq 1 then
      o := Order(chi);
      if o in orders_to_check and not IsDefined(evens_by_order, o) then
        evens_by_order[o] := chi;
      end if;
    end if;
  end for;
  for o -> chi in evens_by_order do
    if o eq 1 then continue; end if;
    d := Dimension(CuspidalSubspace(ModularForms(chi, 2)));
    if d gt best then best := d; end if;
  end for;
  return best;
end function;

function PrecisionTierFor(needed, MIN_PRECISION, MAX_PRECISION)
  tier := MIN_PRECISION;
  while tier lt needed and tier lt MAX_PRECISION do
    tier *:= 2;
  end while;
  if tier gt MAX_PRECISION then tier := MAX_PRECISION; end if;
  return tier;
end function;

function NeededPrecision(dim_estimate, MIN_PRECISION, MAX_PRECISION)
  target_coeffs := Ceiling(MARGIN_MULT_NUM * dim_estimate / MARGIN_MULT_DEN) + BASE_MARGIN_ADD;
  raw := Ceiling(PRECISION_SCALE_NUM * target_coeffs^2 / PRECISION_SCALE_DEN);
  return PrecisionTierFor(Max(raw, MIN_PRECISION), MIN_PRECISION, MAX_PRECISION);
end function;

// ================================================================
// per-ideal processing (unchanged from before, just factored into a
// procedure so it can be called once per tier's shared graded ring)
// ================================================================

procedure ProcessIdeal(N, level_label, chis, GRing, bbs, field_str, ~classical_cache)
  for chi_idx -> chi in chis do
    t0 := Cputime();

    // fields common to every row emitted for this (N, chi_idx)
    base := BlankRow();
    base[C_FIELD] := S(field_str);
    base[C_LEVEL_LABEL] := S(level_label);
    base[C_LEVEL_NORM] := S(Norm(N));
    base[C_CHI_IDX] := S(chi_idx);
    base[C_NUM_CHIS] := S(#chis);

    try
      Mk := HMFSpace(GRing, N, k, chi);
      Dk := DihedralBasis(Mk);
    catch e
      row := base;
      row[C_STATUS] := "err_dihedral_basis";
      row[C_ELAPSED] := S(Cputime()-t0);
      row[C_ERROR] := S(e);
      WriteRow(row);
      continue;
    end try;
    if #Dk eq 0 then continue; end if;

    chi_ord := Order(chi);
    chi_lab := HeckeCharLabel(chi : full_label:=false);
    base[C_CHI_ORDER] := S(chi_ord);
    base[C_CHI_LABEL] := S(chi_lab);
    base[C_NUM_F] := S(#Dk);

    for f_idx -> f in Dk do
      base[C_F_IDX] := S(f_idx);
      for bb_idx -> bb in bbs do
        t0 := Cputime();
        row := base;
        row[C_BB_IDX] := S(bb_idx);

        try
          coeffs := RestrictionToDiagonal(f, GRing, bb : AsCoefficients := true);
        catch e
          row[C_STATUS] := "err_restriction";
          row[C_ELAPSED] := S(Cputime()-t0);
          row[C_ERROR] := S(e);
          WriteRow(row);
          continue;
        end try;

        NN := Level(f);
        classical_level := Integers()!(Denominator(NN)^(-1)*Generator((Denominator(NN)*NN) meet Integers()));
        row[C_CLASSICAL_LEVEL] := S(classical_level);
        row[C_NUM_COEFFS] := S(#coeffs);

        if forall{c : c in coeffs | c eq 0} then
          row[C_STATUS] := "zero_restriction";
          row[C_ELAPSED] := S(Cputime()-t0);
          WriteRow(row);
          continue;
        end if;

        is_rat := true;
        for c in coeffs do
          if not IsCoercible(RationalField(), c) then is_rat := false; break; end if;
        end for;
        if not is_rat then
          row[C_COEFFS_INT] := "false";
          row[C_STATUS] := "non_rational_coeffs";
          row[C_ELAPSED] := S(Cputime()-t0);
          WriteRow(row);
          continue;
        end if;

        target := [RationalField()!c : c in coeffs];
        is_int := forall{c : c in target | Denominator(c) eq 1};
        row[C_COEFFS_INT] := S(is_int);
        NMAX := #target - 1;
        v := [target[idx+1] : idx in [1..NMAX]];

        try
          TryIdentify(classical_level, chi_ord, v, NMAX, ~classical_cache,
            ~status, ~neb_order, ~neb_char_idx, ~classical_dim, ~num_orbits,
            ~matched_orbit_idx, ~minpoly_K, ~deg_K, ~minpoly_c, ~contrib, ~deg_c, ~trace_c, ~num_aut_K);
        catch e
          row[C_STATUS] := "err_identify";
          row[C_ELAPSED] := S(Cputime()-t0);
          row[C_ERROR] := S(e);
          WriteRow(row);
          continue;
        end try;

        row[C_STATUS] := S(status);
        row[C_NEB_ORDER] := S(neb_order);
        row[C_NEB_CHAR_IDX] := S(neb_char_idx);
        row[C_CLASSICAL_DIM] := S(classical_dim);
        row[C_NUM_ORBITS] := S(num_orbits);
        row[C_MATCHED_ORBIT_IDX] := S(matched_orbit_idx);
        row[C_DEG_K] := S(deg_K);
        row[C_MINPOLY_K] := S(minpoly_K);
        row[C_NUM_AUT_K] := S(num_aut_K);
        row[C_MINPOLY_C] := S(minpoly_c);
        row[C_DEG_C] := S(deg_c);
        row[C_TRACE_C] := S(trace_c);

        // best-effort (M, psi) match -- only meaningful if we got a clean single-orbit match
        try
          psi_found, psi_ord, psi_cond_norm, M_poly, M_disc, M_tc, M_classnum := TryMatchGrossenchar(Mk, f);
        catch e
          psi_found := false; psi_ord := 0; psi_cond_norm := 0; M_poly := ""; M_disc := 0;
          M_tc := false; M_classnum := 0;
        end try;
        row[C_PSI_FOUND] := S(psi_found);
        row[C_PSI_ORDER] := S(psi_ord);
        row[C_PSI_COND_NORM] := S(psi_cond_norm);
        row[C_M_POLY] := S(M_poly);
        row[C_M_DISC] := S(M_disc);
        row[C_M_TC] := S(M_tc);
        row[C_M_CLASSNUM] := S(M_classnum);

        if status eq "matched" then
          GetClassicalSpace(classical_level, neb_order, neb_char_idx, ~classical_cache, ~nf_lookup, ~dim_lookup);
          g := nf_lookup[matched_orbit_idx][1];
          Kg := Parent(Coefficient(g,1));
          small_primes := [p : p in PrimesUpTo(40) | classical_level mod p ne 0];
          if Kg cmpeq Rationals() then
            traces := [Coefficient(g,p) : p in small_primes];
          else
            traces := [Trace(Coefficient(g,p)) : p in small_primes];
          end if;
          row[C_HECKE_TRACES] := SerList(traces);
          row[C_CANDIDATE_LABEL] := Sprintf("%o.2.o%oi%o.orb%o", classical_level, neb_order, neb_char_idx, matched_orbit_idx);
        end if;

        row[C_ELAPSED] := S(Cputime()-t0);
        WriteRow(row);
      end for;
    end for;
  end for;
end procedure;

// ================================================================
// main sweep -- phase 1: cheap candidate collection + tier assignment
// (no graded ring touched yet)
// ================================================================

Log(Sprintf("=== starting survey: field=Qsqrt%o, norm range [%o,%o] ===", D, NORM_MIN, LEVEL_BOUND));
if ADAPTIVE_PRECISION then
  Log(Sprintf("adaptive precision: MIN=%o MAX=%o", MIN_PRECISION, MAX_PRECISION));
else
  Log(Sprintf("fixed precision: %o", FIXED_PRECISION));
end if;
t_start := Cputime();

field_str := "Qsqrt" cat IntegerToString(D);
ideals := [N : N in IdealsUpTo(LEVEL_BOUND, F) | not IsZero(N) and Norm(N) ge NORM_MIN];
Log(Sprintf("enumerated %o ideals", #ideals));

// tier (integer precision) -> SeqEnum of <N, level_label, chis>
tiers := AssociativeArray();
num_skipped_no_chis := 0;

for N in ideals do
  level_label := Join(Split(LMFDBLabel(N), "."), "_");
  H := HeckeCharacterGroup(N, [1,2]);
  chis := [chi : chi in Elements(H) | IsCompatibleWeight(chi, k)];
  if #chis eq 0 then
    num_skipped_no_chis +:= 1;
    continue;
  end if;

  if ADAPTIVE_PRECISION then
    chi_orders := {Order(chi) : chi in chis};
    Lest := Norm(N); // classical_level always divides Norm(N) in every case we've observed
    try
      dim_est := EstimateClassicalDim(Lest, chi_orders);
      tier := NeededPrecision(dim_est, MIN_PRECISION, MAX_PRECISION);
    catch e
      tier := MAX_PRECISION; // dimension lookup itself failed -- be conservative
    end try;
  else
    tier := FIXED_PRECISION;
  end if;

  if not IsDefined(tiers, tier) then
    tiers[tier] := [* *]; // List, not SeqEnum -- holds heterogeneous <ideal, string, seq> tuples
  end if;
  Append(~tiers[tier], <N, level_label, chis>);
end for;

Log(Sprintf("phase 1 done (%o s): %o ideals with >=1 compatible chi across %o precision tier(s), %o skipped (no compatible chi)",
  Cputime()-t_start, #ideals - num_skipped_no_chis, #Keys(tiers), num_skipped_no_chis));
for tier -> group in tiers do
  Log(Sprintf("  tier %o: %o ideals", tier, #group));
end for;

// ================================================================
// main sweep -- phase 2: process each tier against its own graded ring
// ================================================================

sorted_tiers := Sort(SetToSequence(Keys(tiers)));
for tier in sorted_tiers do
  group := tiers[tier];
  Log(Sprintf("=== starting tier %o (%o ideals), elapsed so far %o s ===", tier, #group, Cputime()-t_start));
  GRing := GradedRingOfHMFs(F, tier);
  bbs := NarrowClassGroupReps(GRing);
  for entry in group do
    N, level_label, chis := Explode(entry);
    ProcessIdeal(N, level_label, chis, GRing, bbs, field_str, ~classical_cache);
    Log(Sprintf("done level %o (norm %o, tier %o), elapsed total %o s", level_label, Norm(N), tier, Cputime()-t_start));
  end for;
end for;

Log(Sprintf("=== survey complete, total elapsed %o s ===", Cputime()-t_start));
quit;
