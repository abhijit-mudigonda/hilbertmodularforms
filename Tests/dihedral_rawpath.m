/***************************************************
* Regression tests for the raw-GrpHeckeElt weight-[1,1] fast path
* (RawDihedralCandidates / TrivialDiagonalRestrictionGrossenchars /
* PrunedGrossencharsSetRaw / MaximalIdealsOfNormDividing), as distinct
* from Tests/dihedral.m and Tests/_dihedral_basis.m, which exercise the
* older PossibleGrossencharsOfRelQuadExt/DihedralBasis pipeline and never
* touch this code path at all.
*
* Candidate counts below were validated by direct comparison against a
* snapshot of the pipeline taken before the per-divisor modulus
* restructuring (i.e. against the single-large-modulus search, which
* itself inherited only the earlier SetEnum/ConjugateIdeal perf fixes and
* the RM/mixed-signature parity fix, not any of the three correctness
* bugs the restructuring introduced and fixed) -- not just candidate
* counts but, for the two most heavily scrutinized cases, the exact
* multiset of (D_norm, mm_norm, chi_order) too.
***************************************************/

// (m, Ntgt, expected #RawDihedralCandidates) -- spans CM, RM, and
// mixed-signature auxiliary fields; several of these were specifically
// how the trivial-character-duplication, Modulus()-truncation, and
// under-powered-self-conjugacy-test bugs were originally caught.
counts := [
  <2, 15, 3>, <2, 20, 3>, <2, 24, 3>, <2, 30, 9>, <2, 38, 6>, <2, 40, 11>, <2, 48, 12>,
  <3, 26, 9>, <3, 32, 12>, <3, 45, 17>,
  <29, 20, 14>, <29, 21, 4>, <29, 35, 36>,
  <33, 14, 8>, <33, 21, 17>, <33, 26, 44>, <33, 37, 44>,
  <41, 11, 1>
];

for tup in counts do
  m, Ntgt, expected := Explode(tup);
  F := QuadraticField(m);
  ZF := Integers(F);
  GRing := GradedRingOfHMFs(F, 100);
  frakN := ideal<ZF | Ntgt>;
  Dk := RawDihedralCandidates(GRing, frakN, Ntgt);
  assert #Dk eq expected;
end for;

/***************************************************
* Q(sqrt6), N=49: the flagship case that exposed all three modulus-
* restructuring bugs at once (a self-conjugate trivial character
* divides several searched ideals simultaneously; a non-self-conjugate
* pair of primitive characters share a finite conductor but have
* different N_oo; the RM survivor is only correctly recognized as
* non-self-conjugate when tested against the full N_oo of K_abs, not its
* own primitivization-shrunken N_oo). Checked against the exact
* (D_norm, mm_norm, K_totally_complex, chi_order) multiset, not just the
* total count.
***************************************************/

F := QuadraticField(6);
ZF := Integers(F);
GRing := GradedRingOfHMFs(F, 600);
frakN := ideal<ZF | 49>;
Dk := RawDihedralCandidates(GRing, frakN, 49);
assert #Dk eq 15;

rows := [];
for tup in Dk do
  f, D, K, mm := Explode(tup);
  Append(~rows, <Norm(D), Norm(mm), IsTotallyComplex(AbsoluteField(K)), Order(Character(Parent(f)))>);
end for;
rows := Sort(rows);

expected_rows := Sort([
  <2401, 1, false, 2>,
  <2401, 1, true, 2>, <2401, 1, true, 2>, <2401, 1, true, 2>,
  <2401, 1, true, 2>, <2401, 1, true, 2>, <2401, 1, true, 2>,
  <2401, 1, true, 2>, <2401, 1, true, 2>,
  <2401, 1, true, 14>, <2401, 1, true, 14>, <2401, 1, true, 14>,
  <2401, 1, true, 14>, <2401, 1, true, 14>, <2401, 1, true, 14>
]);
assert rows eq expected_rows;

/***************************************************
* RM/CM/mixed-signature spot checks: confirms the weight-[1,1] parity
* filter (fix(Dihedral): filter RM/mixed-signature psi failing
* weight-[1,1] parity) is still correctly wired into the raw-modulus
* search path.
***************************************************/

// genuinely parity-invalid: every candidate here must be excluded, not crash
F41 := QuadraticField(41);
ZF41 := Integers(F41);
G41 := GradedRingOfHMFs(F41, 100);
assert #RawDihedralCandidates(G41, ideal<ZF41|17>, 17) eq 0;

// RM-sourced candidate survives correctly (not spuriously excluded)
F2 := QuadraticField(2);
ZF2 := Integers(F2);
G2 := GradedRingOfHMFs(F2, 100);
ans2 := RawDihedralCandidates(G2, ideal<ZF2|15>, 15);
assert #ans2 eq 3;
found_rm := false;
for tup in ans2 do
  f, D, K, mm := Explode(tup);
  if Norm(D) eq 225 and not IsTotallyComplex(AbsoluteField(K)) then
    found_rm := true;
  end if;
end for;
assert found_rm;

// mixed-signature auxiliary field (signature (2,1) over F=Q(sqrt5))
F5 := QuadraticField(5);
ZF5 := Integers(F5);
G5 := GradedRingOfHMFs(F5, 100);
assert #RawDihedralCandidates(G5, ideal<ZF5|19>, 19) eq 3;
