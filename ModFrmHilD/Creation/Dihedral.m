import "CuspFormFromEigs.m" : codifferent_generator;

/////////////////////////////// Compute quadratic extensions with conductor

intrinsic QuadraticExtensionsWithConductor(NN::RngOrdIdl, InfinityModulus::SeqEnum[RngIntElt] : Dividing := true)
  -> SeqEnum[FldAlg]
  {
    Naive!  Returns the set of quadratic field extensions of conductor equal to NN*oo where 
    oo is InfinityModulus, indexing a subset of real places (as Magma numbers them)
    as in their class field theory machinery.  Use Dividing to allow those with 
    conductor dividing NN*oo.
  }
  ZZF := Order(NN);
  F := NumberField(ZZF);
  _<x> := PolynomialRing(F);
  oo := InfinityModulus;

  U, m := SUnitGroup(NN);
  Usq, msq := quo<U | [2*u : u in Generators(U)]>;
  Ks := [];
  for usq in Usq do
  if usq eq Id(Usq) then continue; end if;
  delta := m(usq@@msq);
  K := ext<F | x^2-delta>;
  ff, ooff := Conductor(AbelianExtension(K));
  if ff eq NN or (Dividing and IsIntegral(NN/ff) and &and[c in oo : c in ooff]) then
    Append(~Ks, K);
  end if;
  end for;
  return Ks;
end intrinsic;

function QuadraticCharacter(I, K)
  // I::RngOrdIdl - An ideal of a field F
  // K::Fld - A quadratic extension of F
  //
  // Returns the value of the quadratic character evaluated 
  // at I. This is equal to (-1)^(#{inert primes factors of I})
  ZK := Integers(K);
  Fact := Factorization(I); 
  sum_inert := 0;
  for foo in Fact do
    P := foo[1];
    PK := ZK !! P;
    FactPK := Factorization(PK);
    if #FactPK eq 1 and FactPK[1][2] eq 1 then // P is inert in K
      sum_inert := sum_inert + foo[2];
    end if;
  end for;
  return (-1)^(sum_inert);
end function;

//////////////////////////////// Conjugates of Grossencharacters

intrinsic ConjugateIdeal(K::Fld, F::Fld, N::RngOrdIdl) -> RngOrdIdl
  {
    inputs:
      K - An absolute number field
      F - A field such that K/F is a quadratic extension
      N - An ideal of the ring of integers of K
    returns:
      The conjugate of N.
  }
  require IsSubfield(F, K) : "F is not a subfield of K";
  K_rel := RelativeField(F, K);
  require Degree(K_rel) eq 2 : "K is not a quadratic extension of F";
  ZK_rel := Integers(K_rel);
  N_rel := ZK_rel!!N;

  // the nontrivial automorphism of K
  aut := Automorphisms(K_rel)[2];
  conj_gens := [aut(gen) : gen in Generators(N_rel)];
  return Integers(K)!!(ideal<ZK_rel | conj_gens>);
end intrinsic;

intrinsic PrunedGrossencharsSet(X::HMFGrossencharsTorsor, F::Fld : RayClassFilter:=false) -> SetEnum
  {
    input:
      X - A HMFGrossencharTorsor of finite order characters
        over a field K with modulus N.
      F - A number field such that K/F is a quadratic extension
        and N is stable under the Gal(K/F) action
      RayClassFilter - optional, passed through to HMFGrossencharsTorsorSet;
        see its documentation.
    returns:
      The set of characters in chi which are not self-conjugate
      under the K/F action, up to conjugation.
  }
  require IsFiniteOrder(X) : "Cannot compute conjugate pairs of\
      infinite order Grossencharacters";
  K := X`BaseField;
  require IsSubfield(F, K) : "F is not a subfield of K";

  N_f, N_oo := Modulus(X);
  H := HeckeCharacterGroup(N_f, N_oo);

  // If H is trivial then the ray class group is trivial and the trivial
  // character is self-conjugate. Checked against H, not S below, since
  // #S eq 1 can also happen from filtering down to one survivor.
  if #H eq 1 then
    return {};
  end if;

  S := HMFGrossencharsTorsorSet(X : RayClassFilter:=RayClassFilter);

  G, mp := RayClassGroup(N_f, N_oo);
  gens := Generators(G);
  idl_gens := [mp(g) : g in gens];
  lcm_order := LCM([Order(g) : g in gens]); 
  L := CyclotomicField(lcm_order);

  chi_to_evals := AssociativeArray();
  evals_to_chi := AssociativeArray();

  pairs := [];

  for chi in S do
    chi_evals := [StrongCoerce(L, chi(I)) : I in idl_gens]; 
    chi_to_evals[chi`RayClassChar] := chi_evals;
    evals_to_chi[chi_evals] := chi;
  end for;

  chis := S;
  pruned_chis := {};
  while not IsEmpty(chis) do
    chi := Rep(chis);
    chi_evals := chi_to_evals[chi`RayClassChar];
    conj_chi_evals := [StrongCoerce(L, chi(ConjugateIdeal(K, F, I))) : I in idl_gens];
    // include chi in the set of pruned chis if
    // chi isn't self-conjugate
    if chi_evals ne conj_chi_evals then
      Include(~pruned_chis, chi);
    end if;
    // remove chi and its conjugate
    assert evals_to_chi[conj_chi_evals] in chis;
    assert chi in chis;
    Exclude(~chis, evals_to_chi[conj_chi_evals]);
    Exclude(~chis, chi);
  end while;
  return pruned_chis;
end intrinsic;

//////////////////////////////// Computing Grossencharacters

function OrderedPlacesOfCMField(K, F : Precision:=25)
  // K - a number field, CM over F
  // F - a number field
  //
  // returns a SeqEnum[RngIntElt] whose ith entry
  // is the index of the place of F lying underneath
  // the ith infinite place of K (in the orderings given by
  // InfinitePlaces(F) and InfinitePlaces(K), respectively).
  
  // check that K/F is a CM extension
  assert IsTotallyReal(F);
  assert IsTotallyComplex(K);
  assert IsSubfield(F, K);
  assert Degree(K) eq 2*Degree(F);

  n := Degree(F);
  F_places := InfinitePlaces(F);
  K_places := InfinitePlaces(K);
  a := F.1;
  a_embed_dict := AssociativeArray();
  for i in [1 .. n] do
    a_i := Round(10^Precision * Evaluate(a, F_places[i]));
    a_embed_dict[a_i] := i;
  end for;

  lies_over := [];
  for i in [1 .. n] do
    b_i := Round(10^Precision * Evaluate(K!a, K_places[i]));
    Append(~lies_over, a_embed_dict[b_i]);
  end for;
  return lies_over;
end function;

intrinsic HeckeCharWeightFromWeight(K::Fld, F::Fld, k::SeqEnum[RngIntElt]) -> SeqEnum[Tup] 
  {}
  // in parallel weight 1, the infinity type is trivial
  if IsParallel(k) and (k[1] eq 1) then
    r, s := Signature(K);
    return [<0, 0> : _ in [1 .. r+s]];
  else
    k_0 := Max(k);
    lies_over := OrderedPlacesOfCMField(K, F);
    hc_wt := [<(k_0 - k[lies_over[i]]) / 2, (k_0 + k[lies_over[i]] - 2) / 2> : i in [1 .. #k]];
    // if the weight is paritious, all the weight components should be integers
    if IsParitious(k) then
      hc_wt := [<Integers()!tup[1], Integers()!tup[2]> : tup in hc_wt];
    end if;
    return hc_wt;
  end if;
end intrinsic;

intrinsic PossibleGrossencharsOfRelQuadExt(K, N, k_hmf, chi : GRing:=false, AllowImprimitive:=false) -> List
  {
    inputs:
      K - relative quadratic extension with base field F and
        discriminant dividing N
      N - integral ideal of F
      k_hmf - weight of HMFs induced by the desired grossencharacters
      chi - (finite order) Hecke character of F of modulus N
      GRing - optional ModFrmHilDGRng; if given, the (chi-independent) set S
        of candidate Grossencharacters is cached on it, keyed by (N, k_hmf, K)
      AllowImprimitive - optional boolean
    returns:
      Grossencharacters of weight k and conductor N/Disc_(K/F) whose
      restriction to AA_F is chi.

      If the weight is parallel, we remove characters which are invariant
      under conjugation (see the ConjugateIdeal intrinsic) and only
      return one character from each pair of conjugate ideals.

      The grossencharacters returned by this function corresponds to distinct CM
      modular forms after (automorphic) induction.

      If AllowImprimitive, we allow the conductor of the Grossencharacter to 
      strictly divide N/Disc_(K/F). Taking ThetaSeries of the resulting 
      characters gives not only the new induced subspace but also the 
      images under the trivial degeneracy map of old induced subspaces. 
  }
  ZK := Integers(K);
  rel_disc := Discriminant(ZK);

  M := N / rel_disc;
  require IsIntegral(M) : "The discriminant of K/F does not divide the level N";

  M := Integers(AbsoluteField(K))!!M;
  K_abs := AbsoluteField(K);
  k := HeckeCharWeightFromWeight(K_abs, BaseField(K), k_hmf);
  X := cHMFGrossencharsTorsor(K_abs, k, M);

  GF, mF := RayClassGroup(N, [1 .. Degree(BaseField(K))]);

  // S doesn't depend on chi, only on (N, k_hmf, K) -- computing it is the
  // expensive step (O(size of the ray class group of K_abs)), so cache it
  // on GRing when available rather than recomputing per character. In the
  // finite order case, filtering via RayClassFilter is fast enough on its
  // own (see PrunedGrossencharsSet/HMFGrossencharsTorsorSet) that it isn't
  // worth caching S separately per chi, so we skip the GRing cache there.
  if IsFiniteOrder(X) then
    // we define grossencharacters over absolute fields
    // so we need to pass in the base field
    F := BaseField(K);
    // Evaluate(cHMFGrossenchar(X,psi), I) reduces to psi(I)^-1 (raw,
    // native GrpHeckeElt evaluation) when IsFiniteOrder(X), since
    // MarkedCharClassRepEvals and EvaluateNoncompactInfinityType are both
    // trivial in that case (see HMFGrossenchar.m). Solving
    // chi(I)*Norm(I)^(k-1) eq psi(I)^-1*QuadraticCharacter(I,K) for
    // psi(I) gives the target values below.
    ray_class_filter := [<Integers(K_abs)!!mF(g),
        QuadraticCharacter(mF(g), K) * chi(mF(g))^-1 * Norm(mF(g))^-(Max(k_hmf) - 1)>
        : g in Generators(GF)];
    S := PrunedGrossencharsSet(X, F : RayClassFilter:=ray_class_filter);
  else
    use_cache := GRing cmpne false;
    k_key := Sprint(k_hmf);
    if use_cache and IsDefined(GRing`GrossencharsSetCache, N)
        and IsDefined(GRing`GrossencharsSetCache[N], k_key)
        and IsDefined(GRing`GrossencharsSetCache[N][k_key], K) then
      S := GRing`GrossencharsSetCache[N][k_key][K];
    else
      S := HMFGrossencharsTorsorSet(X);
      if use_cache then
        if not IsDefined(GRing`GrossencharsSetCache, N) then
          GRing`GrossencharsSetCache[N] := AssociativeArray();
        end if;
        if not IsDefined(GRing`GrossencharsSetCache[N], k_key) then
          GRing`GrossencharsSetCache[N][k_key] := AssociativeArray();
        end if;
        GRing`GrossencharsSetCache[N][k_key][K] := S;
      end if;
    end if;
  end if;

  ans := [* *];
  for psi in S do
    N_psi := ZK!!(Conductor(psi));
    valid_level := (not AllowImprimitive) select \
                 (Norm(N_psi) * rel_disc eq N) else \
                 (IsIntegral(N / (Norm(N_psi) * rel_disc)));
    if valid_level then
      flag := true;
      for g in Generators(GF) do
        I := mF(g);
        if IsAlgebraic(X) then
          nonpar_hack := false;
          // the default value of gen
          gen := 0;
        else
          nonpar_hack := true;
          // in the nonparitious case, we want to feed Evaluate
          // a generator of I which lies in F
          b, x := IsPrincipal(I);
          assert b;
          gen:=x;
        end if;
        flag and:= StrongEquality(
            chi(I) * Norm(I)^(Max(k_hmf) - 1),
            Evaluate(psi, Integers(K_abs)!!(I) : gen:=gen, nonpar_hack:=nonpar_hack) * QuadraticCharacter(I, K)
            );
      end for;
      // S was already filtered by RayClassFilter above in the finite
      // order case, so this recheck should always agree; keep it as an
      // assert-based safety net rather than dropping it, since S is now
      // small this costs nothing.
      if IsFiniteOrder(X) then
        assert flag;
        Append(~ans, psi);
      elif flag then
        Append(~ans, psi);
      end if;
    end if;
  end for;
  return ans;
end intrinsic;

intrinsic PossibleGrossenchars(Mk::ModFrmHilD : AllowImprimitive:=false) -> List
  {
    Given a space Mk of HMFs, computes the grossencharacters which induce
    forms in Mk.

    The induced forms will span the dihedral forms in Mk if Weight(Mk) is parallel
    and the CM forms otherwise.

    If AllowImprimitive, also returns grossencharacters whose conductor
    strictly divides Level(Mk); see PossibleGrossencharsOfRelQuadExt.
  }
  ans := [* *];
  N := Level(Mk);
  GRing := Parent(Mk);
  if IsDefined(GRing`QuadExtWithConductorCache, N) then
    Ks := GRing`QuadExtWithConductorCache[N];
  else
    Ks := QuadraticExtensionsWithConductor(N, [1 .. Degree(BaseField(Mk))]);
    GRing`QuadExtWithConductorCache[N] := Ks;
  end if;

  if not (IsParallel(Weight(Mk)) and Weight(Mk)[1] eq 1) then
    Ks := [K : K in Ks | IsTotallyComplex(K)];
  end if;

  k := Weight(Mk);
  chi := Character(Mk);
  for K in Ks do
    ans cat:= [* psi : psi in PossibleGrossencharsOfRelQuadExt(K, N, k, chi : GRing := GRing, AllowImprimitive := AllowImprimitive) *];
  end for;
  return ans;
end intrinsic;

//////////////////////////////// Computing spaces of dihedral forms

intrinsic MaximalIdealsOfNormDividing(ZK::RngOrd, Mideal::RngOrdIdl) -> SeqEnum
  {
    ZK the ring of integers of a quadratic extension K of F := NumberField
    (Order(Mideal)). Returns a duplicate-free sequence of ideals C of ZK,
    each maximal (under divisibility) among ideals whose relative norm
    N_KF(C) divides Mideal, such that EVERY ideal D of ZK with
    N_KF(D) | Mideal divides some C in the returned sequence (hence
    appears, possibly imprimitively, in HeckeCharacterGroup(C,...) for
    that C).

    This is a much tighter search bound than the single ideal Mideal*ZK
    used previously: since D | N_KF(D)*ZK always (as D*D^sigma =
    N_KF(D)*ZK), Mideal*ZK is a valid but wasteful upper bound whose
    OWN relative norm is Mideal^2, not Mideal -- e.g. for a case with
    #HeckeCharacterGroup(Mideal*ZK,...) = 115248, the ideals actually
    produced by this intrinsic have combined character-group size 146,
    a ~790x reduction, empirically confirmed to translate directly into
    proportionally faster downstream search.

    Per-prime local structure (p a prime of F dividing Mideal with
    exponent e): if p is INERT in K (single prime P, N(P)=p^2), the
    maximal local ideal is P^(e div 2) -- note e need not be even; if e
    is odd, this maximal ideal has norm p^(e-1), strictly less than
    p^e, which is correct (norm p^e is not achievable by any ideal above
    an inert prime at all, so nothing is lost by not trying to hit it
    exactly). If p is RAMIFIED (P^2=p*ZK, N(P)=p), the unique maximal
    ideal is P^e (norm exactly p^e). If p is SPLIT (P,P' distinct,
    N(P)=N(P')=p), every P^a*P'^(e-a) for a in [0..e] is simultaneously
    maximal (norm exactly p^e, but not comparable to each other under
    divisibility) -- e+1 genuinely different ideals to search separately.

    Returns [1*ZK] when Mideal is the unit ideal (vacuously correct: the
    only ideal of norm dividing 1*ZF is 1*ZK itself).
  }
  ZF := Order(Mideal);
  local_choices := [* *];
  for pe in Factorization(Mideal) do
    p := pe[1]; e := pe[2];
    pK_fact := Factorization(ideal<ZK | [ZK!x : x in Generators(p)]>);
    options := [];
    if #pK_fact eq 1 and pK_fact[1][2] eq 2 then
      // ramified: unique maximal ideal, norm exactly p^e
      P := pK_fact[1][1];
      Append(~options, P^e);
    elif #pK_fact eq 1 and pK_fact[1][2] eq 1 then
      // inert: unique maximal ideal, norm p^(2*(e div 2)) <= p^e
      P := pK_fact[1][1];
      Append(~options, P^(e div 2));
    elif #pK_fact eq 2 then
      // split: e+1 pairwise-incomparable maximal ideals, each norm exactly p^e
      P := pK_fact[1][1]; Pp := pK_fact[2][1];
      for a in [0 .. e] do
        Append(~options, P^a * Pp^(e-a));
      end for;
    end if;
    Append(~local_choices, options);
  end for;

  Cs := [1*ZK];
  for opts in local_choices do
    newCs := [];
    for c in Cs do
      for o in opts do
        Append(~newCs, c*o);
      end for;
    end for;
    Cs := newCs;
  end for;
  return Cs;
end intrinsic;

intrinsic TrivialDiagonalRestrictionGrossenchars(K::Fld, frakN::RngOrdIdl, Ntgt::RngIntElt) -> List
  {
    Weight [1,1] fast path: K a quadratic extension of the totally real
    field F with discriminant dividing frakN, Ntgt the target classical
    (rational) level. Returns the (duplicate-free) sequence of finite-order
    psi::GrpHeckeElt whose induced nebentypus is trivial on (Z/NtgtZ)^*,
    i.e. whose diagonal restriction has trivial classical nebentypus.
    Native GrpHeckeElt evaluation throughout; never builds an
    HMFGrossencharsTorsor.

    Searches HeckeCharacterGroup(C,...) across the (few) ideals C returned
    by MaximalIdealsOfNormDividing(ZK, frakN/Discriminant(K)), rather than
    a single HeckeCharacterGroup at the much larger pushforward ideal
    (frakN/Discriminant(K))*ZK -- see that intrinsic's docstring for why
    this is both correct (nothing is missed) and much faster (empirically
    ~790x fewer characters to touch in one measured case).

    Returns a List (not a SeqEnum, since different C give psi from
    different, incompatible HeckeCharacterGroup parents -- a SeqEnum
    requires a common universe) and not a SetEnum (GrpHeckeElt has no
    efficient Hash, so a SetEnum of them is O(n) per Include -- O(n^2)
    total). Every C in Cs is searched (no attempt is made here to skip
    "redundant" conjugate C's -- an earlier version tried that, on the
    reasoning that the filter condition is conjugation-invariant so C and
    C^sigma contribute the same characters up to conjugation, but this is
    unsound: a psi discovered via C need not have ITS OWN true conductor
    equal to C, and that true conductor's conjugate can fail to divide
    C^sigma at all while still dividing some OTHER searched C' -- so
    skipping C^sigma outright can silently drop psi whose conjugate
    partner was never a raw survivor of anything searched. Simpler and
    safer to search everything here and let PrunedGrossencharsSetRaw,
    which handles arbitrary redundancy/overlap correctly, do all the
    conjugate-pairing and self-conjugate exclusion). The resulting raw
    list is not yet deduplicated in any sense -- in particular the SAME
    primitive character can appear more than once, once per C its true
    conductor divides (e.g. the trivial character, conductor 1*ZK,
    divides every C) -- see PrunedGrossencharsSetRaw for how this, and
    self-conjugate exclusion, are handled downstream.
  }
  F := BaseField(K);
  ZF := Integers(F);
  ZK := Integers(K);
  rel_disc := Discriminant(ZK);
  Mideal := frakN / rel_disc;
  require IsIntegral(Mideal) : "The discriminant of K/F does not divide frakN";
  K_abs := AbsoluteField(K);
  N_oo := [1 .. #RealPlaces(K_abs)];

  Uz, mUz := UnitGroup(Integers(Ntgt));
  a_gens := [Integers()!(g @ mUz) : g in Generators(Uz)];
  filter := [<Integers(K_abs)!!ideal<ZF|a>, QuadraticCharacter(ideal<ZF|a>, K)> : a in a_gens];

  Cs := MaximalIdealsOfNormDividing(ZK, Mideal);
  Cs_abs := [Integers(K_abs) !! C : C in Cs];

  S := [* *];
  for C_abs in Cs_abs do
    H := HeckeCharacterGroup(C_abs, N_oo);
    for psi in Elements(H) do
      if &and[psi(t[1]) eq t[2] : t in filter] then
        Append(~S, psi);
      end if;
    end for;
  end for;
  return S;
end intrinsic;

intrinsic RawDihedralCandidates(GRing::ModFrmHilDGRng, frakN::RngOrdIdl, Ntgt::RngIntElt) -> List
  {
    Weight [1,1] fast path replacing DihedralBasis/PossibleGrossenchars/
    HMFGrossenchar for the specific case of searching for forms whose
    diagonal restriction has trivial classical nebentypus at level Ntgt.

    For each auxiliary quadratic extension K of QuadraticExtensionsWithConductor
    (frakN,...), finds every finite-order psi whose induced nebentypus is
    trivial on (Z/NtgtZ)^*, dedupes conjugate pairs, and realizes each as
    an HMF at its own native conductor D | frakN via a raw ThetaSeries,
    lifted to frakN via Inclusion when D properly divides frakN.

    Returns a List of tuples <f, D, K, mm> where f is the level-frakN HMF,
    D is psi's native conductor, K is the auxiliary field, mm := frakN/D
    is the degeneracy scalar.
  }
  F := BaseField(GRing);
  ZF := Integers(F);
  Ks := QuadraticExtensionsWithConductor(frakN, [1 .. Degree(F)]);

  ans := [* *];
  for K in Ks do
    S := TrivialDiagonalRestrictionGrossenchars(K, frakN, Ntgt);
    if #S eq 0 then continue; end if;
    K_abs := AbsoluteField(K);
    S := PrunedGrossencharsSetRaw(S, K_abs, F);

    ZK := Integers(K);
    rel_disc := Discriminant(ZK);
    for psi0 in S do
      psi := AssociatedPrimitiveCharacter(psi0);
      D := ZF !! (Norm(ZK!!Conductor(psi)) * rel_disc);
      // psi's native conductor may not divide frakN at all (e.g. a
      // psi primitive at the full search modulus M=frakN/rel_disc has
      // D = Norm(M*O_K)*rel_disc = frakN^2/rel_disc, which generally
      // doesn't divide frakN) -- such psi are not realizable at this
      // level, old or new, and must be skipped rather than forced
      // through Inclusion.
      if not IsIntegral(frakN/D) then
        continue;
      end if;

      GD, mpD := RayClassGroup(D, [1 .. Degree(F)]);
      B_D := <<mpD(g), psi(Integers(K_abs)!!mpD(g))^-1 * QuadraticCharacter(mpD(g), K)> : g in Generators(GD)>;
      chi_D := HeckeCharacter(D, [1 .. Degree(F)], B_D);
      // chi_D is reconstructed directly from psi's values, with no guarantee
      // it's an admissible weight-[1,1] nebentype: unlike the old pipeline
      // (which only ever tries chi's that are already elements of a
      // HeckeCharacterGroup, hence automatically compatible), we build chi_D
      // from scratch and psi's real-place sign behavior can make it
      // incompatible -- concretely, at a real place v of F that splits into
      // two real places w,w' of K (i.e. K is unramified/real at v: this is
      // the generic case in the RM branch, and can also occur at individual
      // places of a mixed-signature K), weight-[1,1] compatibility requires
      // psi's local signs at w and w' to differ; nothing upstream enforces
      // this. Skip (not error) when it fails, matching what the old pipeline
      // would implicitly do (no compatible chi ever gets tried against this
      // psi, so it's silently absent from old's results too).
      if not IsCompatibleWeight(chi_D, [1,1]) then
        continue;
      end if;
      Mk_D := HMFSpace(GRing, D, [1,1], chi_D);
      f_D := ThetaSeries(Mk_D, psi);

      if D ne frakN then
        GN, mpN := RayClassGroup(frakN, [1 .. Degree(F)]);
        B_N := <<mpN(g), psi(Integers(K_abs)!!mpN(g))^-1 * QuadraticCharacter(mpN(g), K)> : g in Generators(GN)>;
        chi_N := HeckeCharacter(frakN, [1 .. Degree(F)], B_N);
        // chi_N is built from the same psi via the same construction, so
        // expected to inherit compatibility from chi_D, but that's not
        // structurally guaranteed (different modulus, different generating
        // set) -- check explicitly rather than assume.
        if not IsCompatibleWeight(chi_N, [1,1]) then
          continue;
        end if;
        Mk_N := HMFSpace(GRing, frakN, [1,1], chi_N);
        // Inclusion(f_D,Mk_N) (2-arg) returns one form per divisor dd of
        // frakN/D, each computed via the 3-arg Inclusion(f_D,Mk_N,dd) using
        // dd (not the full frakN/D) as its own individual degeneracy scale.
        // Loop over the divisors ourselves so each returned form gets
        // tagged with the dd actually used to build it, not the fixed
        // frakN/D -- tagging every form with frakN/D regardless of which
        // specific divisor produced it is wrong provenance metadata.
        for dd in Divisors(frakN/D) do
          f := Inclusion(f_D, Mk_N, dd);
          Append(~ans, <f, D, K, dd>);
        end for;
      else
        Append(~ans, <f_D, D, K, 1*ZF>);
      end if;
    end for;
  end for;
  return ans;
end intrinsic;

intrinsic DebugRawDihedralCandidatesInstrumented(GRing::ModFrmHilDGRng, frakN::RngOrdIdl, Ntgt::RngIntElt) -> List, Assoc
  {
    TEMPORARY debug utility: byte-for-byte copy of RawDihedralCandidates
    with Cputime() checkpoints inserted, to find discrepancies between
    isolated component timing and the real end-to-end cost (everything
    computed fresh in one call, no reused objects across separate script
    loads that might benefit from stale caching).
  }
  F := BaseField(GRing);
  ZF := Integers(F);

  t0 := Cputime();
  Ks := QuadraticExtensionsWithConductor(frakN, [1 .. Degree(F)]);
  timers := AssociativeArray();
  timers["QuadraticExtensionsWithConductor"] := Cputime(t0);
  timers["Trivial"] := 0.0;
  timers["Pruned"] := 0.0;
  timers["AbsoluteField"] := 0.0;
  timers["PerPsiLoop"] := 0.0;
  timers["RCG_D"] := 0.0;
  timers["HC_D"] := 0.0;
  timers["Compat_D"] := 0.0;
  timers["HMFSpace_D"] := 0.0;
  timers["ThetaSeries"] := 0.0;
  timers["RCG_N"] := 0.0;
  timers["HC_N"] := 0.0;
  timers["Compat_N"] := 0.0;
  timers["HMFSpace_N"] := 0.0;
  timers["Inclusion"] := 0.0;
  n_theta := 0;

  ans := [* *];
  for K in Ks do
    t0 := Cputime();
    S := TrivialDiagonalRestrictionGrossenchars(K, frakN, Ntgt);
    timers["Trivial"] +:= Cputime(t0);
    if #S eq 0 then continue; end if;

    t0 := Cputime();
    K_abs := AbsoluteField(K);
    timers["AbsoluteField"] +:= Cputime(t0);

    t0 := Cputime();
    S := PrunedGrossencharsSetRaw(S, K_abs, F);
    timers["Pruned"] +:= Cputime(t0);

    ZK := Integers(K);
    rel_disc := Discriminant(ZK);
    for psi0 in S do
      t_loop0 := Cputime();
      psi := AssociatedPrimitiveCharacter(psi0);
      D := ZF !! (Norm(ZK!!Conductor(psi)) * rel_disc);
      if not IsIntegral(frakN/D) then
        timers["PerPsiLoop"] +:= Cputime(t_loop0);
        continue;
      end if;

      t0 := Cputime();
      GD, mpD := RayClassGroup(D, [1 .. Degree(F)]);
      timers["RCG_D"] +:= Cputime(t0);

      t0 := Cputime();
      B_D := <<mpD(g), psi(Integers(K_abs)!!mpD(g))^-1 * QuadraticCharacter(mpD(g), K)> : g in Generators(GD)>;
      chi_D := HeckeCharacter(D, [1 .. Degree(F)], B_D);
      timers["HC_D"] +:= Cputime(t0);

      t0 := Cputime();
      is_compat := IsCompatibleWeight(chi_D, [1,1]);
      timers["Compat_D"] +:= Cputime(t0);
      if not is_compat then
        timers["PerPsiLoop"] +:= Cputime(t_loop0);
        continue;
      end if;

      t0 := Cputime();
      Mk_D := HMFSpace(GRing, D, [1,1], chi_D);
      timers["HMFSpace_D"] +:= Cputime(t0);

      t0 := Cputime();
      f_D := ThetaSeries(Mk_D, psi);
      timers["ThetaSeries"] +:= Cputime(t0);
      n_theta +:= 1;

      if D ne frakN then
        t0 := Cputime();
        GN, mpN := RayClassGroup(frakN, [1 .. Degree(F)]);
        timers["RCG_N"] +:= Cputime(t0);

        t0 := Cputime();
        B_N := <<mpN(g), psi(Integers(K_abs)!!mpN(g))^-1 * QuadraticCharacter(mpN(g), K)> : g in Generators(GN)>;
        chi_N := HeckeCharacter(frakN, [1 .. Degree(F)], B_N);
        timers["HC_N"] +:= Cputime(t0);

        t0 := Cputime();
        n_compat_n := IsCompatibleWeight(chi_N, [1,1]);
        timers["Compat_N"] +:= Cputime(t0);
        if not n_compat_n then
          timers["PerPsiLoop"] +:= Cputime(t_loop0);
          continue;
        end if;

        t0 := Cputime();
        Mk_N := HMFSpace(GRing, frakN, [1,1], chi_N);
        timers["HMFSpace_N"] +:= Cputime(t0);

        t0 := Cputime();
        for dd in Divisors(frakN/D) do
          f := Inclusion(f_D, Mk_N, dd);
          Append(~ans, <f, D, K, dd>);
        end for;
        timers["Inclusion"] +:= Cputime(t0);
      else
        Append(~ans, <f_D, D, K, 1*ZF>);
      end if;
      timers["PerPsiLoop"] +:= Cputime(t_loop0);
    end for;
  end for;
  timers["n_theta"] := n_theta;
  timers["n_final"] := #ans;
  return ans, timers;
end intrinsic;

intrinsic DebugPerCandidateTiming(GRing::ModFrmHilDGRng, K::Fld, frakN::RngOrdIdl, Ntgt::RngIntElt) -> Tup
  {
    TEMPORARY debug utility: like the per-K inner loop of
    RawDihedralCandidates, but returns cumulative timings for each stage
    instead of the candidates themselves.
  }
  F := BaseField(GRing);
  ZF := Integers(F);
  S := TrivialDiagonalRestrictionGrossenchars(K, frakN, Ntgt);
  K_abs := AbsoluteField(K);
  S := PrunedGrossencharsSetRaw(S, K_abs, F);
  ZK := Integers(K);
  rel_disc := Discriminant(ZK);

  t_rcg := 0.0; t_hc := 0.0; t_compat := 0.0; t_hmfspace := 0.0; t_theta := 0.0; t_incl := 0.0;
  n_pass := 0; n_compat := 0; n_theta := 0;
  for psi0 in S do
    psi := AssociatedPrimitiveCharacter(psi0);
    D := ZF !! (Norm(ZK!!Conductor(psi)) * rel_disc);
    if not IsIntegral(frakN/D) then continue; end if;
    n_pass +:= 1;

    t0 := Cputime();
    GD, mpD := RayClassGroup(D, [1 .. Degree(F)]);
    t_rcg +:= Cputime(t0);

    t0 := Cputime();
    B_D := <<mpD(g), psi(Integers(K_abs)!!mpD(g))^-1 * QuadraticCharacter(mpD(g), K)> : g in Generators(GD)>;
    chi_D := HeckeCharacter(D, [1 .. Degree(F)], B_D);
    t_hc +:= Cputime(t0);

    t0 := Cputime();
    is_compat := IsCompatibleWeight(chi_D, [1,1]);
    t_compat +:= Cputime(t0);
    if not is_compat then continue; end if;
    n_compat +:= 1;

    t0 := Cputime();
    Mk_D := HMFSpace(GRing, D, [1,1], chi_D);
    t_hmfspace +:= Cputime(t0);

    t0 := Cputime();
    f_D := ThetaSeries(Mk_D, psi);
    t_theta +:= Cputime(t0);
    n_theta +:= 1;

    if D ne frakN then
      t0 := Cputime();
      GN, mpN := RayClassGroup(frakN, [1 .. Degree(F)]);
      B_N := <<mpN(g), psi(Integers(K_abs)!!mpN(g))^-1 * QuadraticCharacter(mpN(g), K)> : g in Generators(GN)>;
      chi_N := HeckeCharacter(frakN, [1 .. Degree(F)], B_N);
      if IsCompatibleWeight(chi_N, [1,1]) then
        Mk_N := HMFSpace(GRing, frakN, [1,1], chi_N);
        for dd in Divisors(frakN/D) do
          f := Inclusion(f_D, Mk_N, dd);
        end for;
      end if;
      t_incl +:= Cputime(t0);
    end if;
  end for;
  return <n_pass, n_compat, n_theta, t_rcg, t_hc, t_compat, t_hmfspace, t_theta, t_incl>;
end intrinsic;

intrinsic RawDihedralPreInclusionForms(GRing::ModFrmHilDGRng, K::Fld, frakN::RngOrdIdl, Ntgt::RngIntElt) -> List
  {
    Debug utility: like RawDihedralCandidates, but restricted to a single
    given K and returning the PRE-Inclusion theta series f_D (at D, its own
    primitive conductor) rather than the post-Inclusion frakN-level lift.
    Returns a List of <f_D, D, chi_D>.
  }
  F := BaseField(GRing);
  ZF := Integers(F);
  S := TrivialDiagonalRestrictionGrossenchars(K, frakN, Ntgt);
  if #S eq 0 then return [* *]; end if;
  K_abs := AbsoluteField(K);
  S := PrunedGrossencharsSetRaw(S, K_abs, F);

  ZK := Integers(K);
  rel_disc := Discriminant(ZK);
  out := [* *];
  for psi0 in S do
    psi := AssociatedPrimitiveCharacter(psi0);
    D := ZF !! (Norm(ZK!!Conductor(psi)) * rel_disc);
    if not IsIntegral(frakN/D) then continue; end if;

    GD, mpD := RayClassGroup(D, [1 .. Degree(F)]);
    B_D := <<mpD(g), psi(Integers(K_abs)!!mpD(g))^-1 * QuadraticCharacter(mpD(g), K)> : g in Generators(GD)>;
    chi_D := HeckeCharacter(D, [1 .. Degree(F)], B_D);
    Mk_D := HMFSpace(GRing, D, [1,1], chi_D);
    f_D := ThetaSeries(Mk_D, psi);
    Append(~out, <f_D, D, chi_D>);
  end for;
  return out;
end intrinsic;

intrinsic RawDihedralPsiInfo(GRing::ModFrmHilDGRng, frakN::RngOrdIdl, Ntgt::RngIntElt) -> List
  {
    Debug utility: like RawDihedralCandidates but returns identifying data
    instead of theta series, for diffing against the old pipeline's
    candidate set. Returns a List of <K_disc_norm, K_tc, D_norm, chi_N>
    where chi_N is the reconstructed nebentypus AT frakN (i.e. after the
    same restriction/lift logic RawDihedralCandidates uses to build Mk_N).
  }
  F := BaseField(GRing);
  ZF := Integers(F);
  Ks := QuadraticExtensionsWithConductor(frakN, [1 .. Degree(F)]);

  out := [* *];
  for K in Ks do
    S := TrivialDiagonalRestrictionGrossenchars(K, frakN, Ntgt);
    if #S eq 0 then continue; end if;
    K_abs := AbsoluteField(K);
    S := PrunedGrossencharsSetRaw(S, K_abs, F);

    ZK := Integers(K);
    rel_disc := Discriminant(ZK);
    for psi0 in S do
      psi := AssociatedPrimitiveCharacter(psi0);
      D := ZF !! (Norm(ZK!!Conductor(psi)) * rel_disc);
      if not IsIntegral(frakN/D) then continue; end if;

      GN, mpN := RayClassGroup(frakN, [1 .. Degree(F)]);
      B_N := <<mpN(g), psi(Integers(K_abs)!!mpN(g))^-1 * QuadraticCharacter(mpN(g), K)> : g in Generators(GN)>;
      chi_N := HeckeCharacter(frakN, [1 .. Degree(F)], B_N);
      Append(~out, <Norm(rel_disc), IsTotallyComplex(K), Norm(D), chi_N>);
    end for;
  end for;
  return out;
end intrinsic;

intrinsic PrunedGrossencharsSetRaw(S::List, K_abs::FldNum, F::Fld) -> List
  {
    Raw-GrpHeckeElt analog of PrunedGrossencharsSet: given a sequence S of
    finite-order psi::GrpHeckeElt of K_abs (drawn from possibly several
    different ambient moduli -- see TrivialDiagonalRestrictionGrossenchars,
    which now searches HeckeCharacterGroup(C,...) across several ideals C
    rather than one large modulus), returns one representative of each
    genuinely distinct, non-self-conjugate primitive character among them,
    i.e. exactly the set that gives distinct cuspidal dihedral theta
    series. Unlike PrunedGrossencharsSet, does not wrap elements in
    HMFGrossenchar and does not use StrongCoerce -- native evaluations of
    elements of the same HeckeCharacterGroup already land in a common
    field.

    Two things S can now contain that the single-shared-modulus version
    this replaces did not need to worry about:

    (1) The SAME primitive character can appear more than once in S, via
    different ambient moduli C that its true conductor happens to divide
    (e.g. the trivial character, conductor 1*ZK, divides every C searched
    and so is a raw survivor of every one of them; more generally any
    character whose conductor divides the "gcd" of two or more searched
    C's). Handled by primitivizing every element up front and deduping on
    the primitive result before any conjugate-pairing logic runs.

    Deduping is keyed on <N_f, N_oo, Eltseq(psi)>, NOT on Modulus(psi)
    used directly as (part of) a tuple: Modulus returns TWO values
    (N_f, N_oo), and embedding a multi-return call directly inside a
    tuple/list-forming expression in Magma silently keeps only the FIRST
    of them -- a real bug caught here: AssociatedPrimitiveCharacter can
    shrink N_oo to a proper subset of the places actually present in the
    ambient modulus (not just the finite part), and two psi with the same
    finite conductor but genuinely different N_oo (e.g. [1,4] vs [2,3],
    swapped by the K/F automorphism) were being wrongly merged as
    "duplicates" when only N_f was captured -- silently discarding one of
    two distinct primitive characters. N_f, N_oo are pulled out via
    explicit assignment first specifically to avoid this.

    (2) A psi that is individually self-conjugate (psi eq psi^sigma --
    these don't correspond to cuspidal dihedral forms, only to
    reducible/Eisenstein-type inductions) can arise even when the ambient
    modulus C it was found via is NOT itself self-conjugate: self-
    conjugacy is a property of the character's own (primitive) conductor,
    not of whichever possibly-non-primitive C it happened to be
    discovered through. So self-conjugacy is checked (and conjugate
    pairing is done) only AFTER primitivizing and deduping per (1) above,
    directly on each element's own true conductor -- never inferred from
    which C it came from.

    The primitivize+dedupe set is small in every case actually
    encountered (bounded by the total size of the -- already small, see
    MaximalIdealsOfNormDividing -- HeckeCharacterGroup(C,...) searched for
    each C), so the O(n^2)-in-the-deduped-set conjugate-pairing loop below
    is not a performance concern; the historical O(n^2) blowup this
    intrinsic was rewritten to avoid (n in the thousands, from a single
    huge shared modulus) cannot recur under the current caller.
  }
  if #S eq 0 then
    return S;
  end if;

  // (1) primitivize + dedupe exact duplicates. Keyed on a string, not a
  // tuple: tuples <Nf, Noo, Eltseq(psi)> vary in "shape" across different
  // psi (Noo and Eltseq(psi) can have different lengths), which breaks
  // SetEnum's homogeneous-universe requirement (same underlying issue as
  // GrpHeckeElt having no efficient Hash -- see TrivialDiagonalRestrictionGrossenchars's
  // docstring); a string is always a uniform, hashable type.
  seen_keys := {};
  prims := [* *];
  for psi0 in S do
    psi := AssociatedPrimitiveCharacter(psi0);
    Nf, Noo := Modulus(psi);
    key := Sprintf("%o|%o|%o", Nf, Noo, Eltseq(psi));
    if key in seen_keys then continue; end if;
    Include(~seen_keys, key);
    Append(~prims, psi);
  end for;

  // (2) exclude self-conjugate psi, pair up the rest. Test generators are
  // drawn from RayClassGroup(Nf_i, N_oo_full) -- the FULL set of infinite
  // places of K_abs, not psi's own (post-primitivization) N_oo -- even
  // though psi eq AssociatedPrimitiveCharacter(psi) may only "genuinely"
  // depend on a proper subset of the places (e.g. N_oo={1,4} out of
  // {1,2,3,4}). This matters: two primitive characters can have the same
  // finite conductor and each look individually self-conjugate when
  // tested only against their own (shrunken) N_oo's ray class group,
  // while genuinely being a conjugate PAIR of each other once tested
  // against the full infinity structure -- confirmed against the old,
  // never-primitivizing-before-comparing pipeline, which evaluates
  // exactly this way (implicitly, by testing while still at the original
  // un-primitivized ambient modulus) and does NOT exclude the analogous
  // survivor. Using a too-small N_oo for the test is a strictly weaker,
  // potentially-wrong test, not just an inefficient one.
  N_oo_full := [1 .. #RealPlaces(K_abs)];
  pruned_chis := [* *];
  excluded := {};
  for i -> psi in prims do
    if i in excluded then continue; end if;
    Nf_i := Modulus(psi);
    Csigma := ConjugateIdeal(K_abs, F, Nf_i);
    G_i, mp_i := RayClassGroup(Nf_i, N_oo_full);
    idl_gens_i := [mp_i(g) : g in Generators(G_i)];
    conj_idl_gens_i := [ConjugateIdeal(K_abs, F, I) : I in idl_gens_i];
    psi_evals := [psi(I) : I in idl_gens_i];

    if Nf_i eq Csigma and psi_evals eq [psi(cI) : cI in conj_idl_gens_i] then
      // self-conjugate: exclude entirely, no partner to find
      Include(~excluded, i);
      continue;
    end if;

    // find the (unique) partner j: primitive, finite conductor Csigma,
    // agreeing with psi under conjugation
    found := false;
    for j -> psi2 in prims do
      if j eq i or j in excluded then continue; end if;
      if Modulus(psi2) ne Csigma then continue; end if;
      if psi_evals eq [psi2(cI) : cI in conj_idl_gens_i] then
        Include(~excluded, j);
        found := true;
        break;
      end if;
    end for;
    assert found;
    Include(~excluded, i);
    Append(~pruned_chis, psi);
  end for;
  return pruned_chis;
end intrinsic;

intrinsic ThetaSeries(Mk::ModFrmHilD, psi::GrpHeckeElt) -> ModFrmHilDElt
  {
    Weight [1,1] fast path: computes the theta series of a finite-order
    Hecke character psi of K directly via native GrpHeckeElt evaluation,
    without wrapping psi in an HMFGrossenchar. psi must already be
    primitive (AssociatedPrimitiveCharacter(psi)) -- Magma's native
    evaluation silently returns 0 off ideals not coprime to psi's
    ambient modulus rather than erroring, so an imprimitive psi would
    silently corrupt coefficients rather than fail loudly.
  }
  require Weight(Mk) eq [1,1] : "ThetaSeries(Mk, psi::GrpHeckeElt) only implements weight [1,1]";
  require Modulus(psi) eq Conductor(psi) : "psi must be primitive; call AssociatedPrimitiveCharacter(psi) first";
  M := Parent(Mk);
  K := NumberField(Order(Modulus(psi)));

  coeffs_by_pp := AssociativeArray();
  for pp in PrimeIdeals(M) do
    fact := Factorization(Integers(K) !! pp);
    g := #fact;
    if g eq 2 then
      coeffs_by_pp[pp] := StrongAdd([* psi(fact[1][1]), psi(fact[2][1]) *]);
    elif fact[1][2] ne 1 then
      coeffs_by_pp[pp] := psi(fact[1][1]);
    else
      coeffs_by_pp[pp] := 0;
    end if;
  end for;

  nonzero_coeff_vals := [* v : pp -> v in coeffs_by_pp | v ne 0 *];
  // FieldOfFractions, not a bare Parent(): if every nonzero coefficient
  // happens to be a plain rational integer, Parent(StrongAdd(...)) infers
  // Integers() (type RngInt), which is a ring but not a Fld -- and
  // CuspEigenformFromCoeffsAtPrimes's coeff_ring downstream (IdlCoeffToEltCoeff)
  // requires a genuine field, crashing with a type mismatch otherwise.
  coeff_ring := (#nonzero_coeff_vals eq 0) select Rationals() else FieldOfFractions(Parent(StrongAdd(nonzero_coeff_vals)));
  return CuspEigenformFromCoeffsAtPrimes(Mk, coeffs_by_pp : coeff_ring := coeff_ring);
end intrinsic;

intrinsic ThetaSeries(Mk::ModFrmHilD, psi::HMFGrossenchar) -> ModFrmHilDElt
  {
    Given a totally real field F, a quadratic extension K of F,
    and a finite order Hecke character of K, compute the associated theta series.
  }
  M := Parent(Mk);
  F := BaseField(M);
 
  ZF := Integers(F);
  prec := Precision(M);
  K := NumberField(Order(Modulus(psi))); 

  // TODO abhijitm this is known (9/20/25) to fail in a lot of nonparitious
  // cases due to field coercion issues. I'm hoping that the pending migration
  // to FldNumEmb will resolve them, but until then use with caution.
  //
  // You might be able to get away with calling #PossibleGrossenchars, which should
  // work more often, and using it as an upper bound for the dimension of CM forms.
  paritious_weight := IsParitious(Weight(Mk));
  if not paritious_weight then
    X := Parent(psi);
    require NarrowClassNumber(F) eq 1 : "Nonparitious forms are only implemented for
      fields with narrow class number one";
    assert IsSubfield(F, K);
    K_rel := RelativeField(F, K);
    // K/F is CM and tau is the nontrivial automorphism
    assert Degree(K_rel) eq 2;
    tau := Automorphisms(K_rel)[2];
    k := Weight(Mk);
    k0 := Max(k);
    lies_over := OrderedPlacesOfCMField(K, F);
    shifted_half_weight := [(k0 - k[lies_over[i]]) / 2 : i in [1..#k]];
    custom_weight := [
      <Integers()!(X`Weight[i][1] - shifted_half_weight[i]),
      Integers()!(X`Weight[i][2] - shifted_half_weight[i])> :
      i in [1 .. Degree(F)]];
    // TODO abhijitm somehow some errors don't occur when this
    // is X`BaseField instead of the splitting field.
    L := SplittingField(X`BaseField);
  end if;

  // We create an associative array indexed by prime ideals pp up to 
  // Precision(Parent(Mk)) and populate them with traces associated to psi.
  coeffs_by_pp := AssociativeArray();
  pis_by_pp := AssociativeArray();
  for pp in PrimeIdeals(M) do
    fact := Factorization(Integers(K) !! pp);
    g := #fact;
    d := InertiaDegree(pp);
    if paritious_weight then
      if g eq 2 then
        coeffs_by_pp[pp] := StrongAdd([* psi(fact[1][1]), psi(fact[2][1]) *]);
      elif fact[1][2] ne 1 then
        coeffs_by_pp[pp] := psi(fact[1][1]);
      else
        coeffs_by_pp[pp] := 0;
      end if;
    else
      if g eq 2 then
        _, lambda_1 := IsNarrowlyPrincipal(fact[1][1]);
        lambdas := [lambda_1, tau(lambda_1)];
        pi := F!(&*lambdas);
        assert Norm(pi) eq Norm(pp);
        psi_evals := [* Evaluate(psi, fact[i][1] : custom_weight:=custom_weight, gen:=lambdas[i])
                  : i in [1..2] *];
        coeffs_by_pp[pp] := L!StrongAdd(psi_evals);
      elif fact[1][2] ne 1 then
        _, lambda := IsNarrowlyPrincipal(fact[1][1]);
        lambdas := [lambda, tau(lambda)];
        pi := F!(&*lambdas);
        assert Norm(pi) eq Norm(pp);
        coeffs_by_pp[pp] := L!Evaluate(psi, fact[1][1] : custom_weight:=custom_weight, gen:=lambda);
      else
        coeffs_by_pp[pp] := L!0;
        pi := codifferent_generator(Parent(Mk))^-1 * IdealToRep(M, pp);
      end if;
      pis_by_pp[pp] := pi;
    end if;
  end for;

  if paritious_weight then
    nonzero_coeff_vals := [* v : pp -> v in coeffs_by_pp | v ne 0 *];
    coeff_ring := (#nonzero_coeff_vals eq 0) select Rationals() else Parent(StrongAdd(nonzero_coeff_vals));
    return CuspEigenformFromCoeffsAtPrimes(Mk, coeffs_by_pp : coeff_ring := coeff_ring);
  else
    return CuspEigenformFromCoeffsAtPrimes(Mk, coeffs_by_pp : 
                                                    from_a_pp:=false,
                                                    mfh_reps:=pis_by_pp,
                                                    coeff_ring:=L
                                                  );
  end if;
end intrinsic;

intrinsic ProbabilisticDihedralTest(f::ModFrmHilDElt) -> BoolElt
  {returns true if this form could be dihedral, false if it cannot be}
  Mk := Parent(f);
  M := Parent(Mk);
  F := BaseField(Mk);
  N := Level(Mk);
  k := Weight(Mk);
  BOUND := 100;
  is_paritious := IsParitious(k);

  Ks := QuadraticExtensionsWithConductor(N, [1 .. Degree(F)]);
  for K in Ks do
    possibly_dihedral := true;
    ZK := Integers(K);

    // inert primes stores the inert primes of norm at most BOUND
    inert_primes := [pp : pp in PrimesUpTo(BOUND, F : coprime_to:=Discriminant(ZK))\
                      | #Factorization(ZK!!pp) eq 1];
    for pp in inert_primes do
      if is_paritious then
        // For paritious forms, use the existing approach
        if not IsZero(Coefficient(f, pp)) then
          possibly_dihedral := false;
          break;
        end if;
      else
        // For nonparitious forms, check if a_nu = 0 where nu is the EltCoeff 
        // corresponding to the ideal pp
        nu := IdealToRep(M, pp);
        bb := IdealToNarrowClassRep(M, pp);
        if not IsZero(Coefficient(f, bb, nu)) then
          possibly_dihedral := false;
          break;
        end if;
      end if;
    end for;
    if possibly_dihedral then
      return true;
    end if;
  end for;
  return false;
end intrinsic;

intrinsic DebugPrunedSetComparison(K::Fld, psi::GrpHeckeElt, k_hmf::SeqEnum[RngIntElt]) -> BoolElt, BoolElt, RngIntElt, RngIntElt, RngOrdIdl, GrpHeckeElt
  {
    Debug utility: given psi already primitive (AssociatedPrimitiveCharacter
    applied), computes its own level D and induced chi_D exactly as
    RawDihedralCandidates/PossibleGrossencharsOfRelQuadExt would, then
    reconstructs the exact internal steps of PossibleGrossencharsOfRelQuadExt's
    finite-order branch for (K,D,chi_D), checking whether psi appears -- via
    its underlying RayClassChar -- in the RayClassFilter-optimized
    PrunedGrossencharsSet result (what old's search actually uses) and/or
    the UNFILTERED PrunedGrossencharsSet result (RayClassFilter bypassed
    entirely). Returns <found_in_filtered, found_in_unfiltered,
    #S_filtered, #S_unfiltered, D>.
  }
  ZK := Integers(K);
  rel_disc := Discriminant(ZK);
  F := BaseField(K);
  ZF := Integers(F);
  K_abs := AbsoluteField(K);
  D := ZF !! (Norm(ZK!!Conductor(psi)) * rel_disc);
  require IsIntegral(D/rel_disc) : "D/rel_disc not integral -- shouldn't happen";

  GD0, mpD0 := RayClassGroup(D, [1 .. Degree(F)]);
  B_D := <<mpD0(g), psi(Integers(K_abs)!!mpD0(g))^-1 * QuadraticCharacter(mpD0(g), K)> : g in Generators(GD0)>;
  chi_D := HeckeCharacter(D, [1 .. Degree(F)], B_D);

  M := D / rel_disc;
  M := Integers(K_abs)!!M;
  k := HeckeCharWeightFromWeight(K_abs, F, k_hmf);
  X := cHMFGrossencharsTorsor(K_abs, k, M);

  GF, mF := RayClassGroup(D, [1 .. Degree(F)]);
  ray_class_filter := [<Integers(K_abs)!!mF(g),
      QuadraticCharacter(mF(g), K) * chi_D(mF(g))^-1 * Norm(mF(g))^-(Max(k_hmf) - 1)>
      : g in Generators(GF)];
  S_filtered := PrunedGrossencharsSet(X, F : RayClassFilter:=ray_class_filter);
  S_unfiltered := PrunedGrossencharsSet(X, F);

  found_filtered := exists{c : c in S_filtered | c`RayClassChar eq psi};
  found_unfiltered := exists{c : c in S_unfiltered | c`RayClassChar eq psi};
  return found_filtered, found_unfiltered, #S_filtered, #S_unfiltered, D, chi_D;
end intrinsic;

intrinsic CompareNewOldDihedralForK(GRing::ModFrmHilDGRng, K::Fld, frakN::RngOrdIdl, Ntgt::RngIntElt) -> List
  {
    Debug/testing utility: for a single auxiliary K, groups the new
    (raw GrpHeckeElt) pipeline's surviving psi by their reconstructed
    chi_N, and for each chi_N compares against the old
    (PossibleGrossencharsOfRelQuadExt/HMFGrossenchar) pipeline restricted
    to the same K. Returns a List of tuples <chi_N, num_new, num_old,
    same_dim, matches> where matches[i] is true iff new_forms[i] equals
    (as a ModFrmHilDElt, via eq) some form in old_forms.
  }
  F := BaseField(GRing);
  ZF := Integers(F);
  ZK := Integers(K);
  rel_disc := Discriminant(ZK);
  K_abs := AbsoluteField(K);

  S := TrivialDiagonalRestrictionGrossenchars(K, frakN, Ntgt);
  S := PrunedGrossencharsSetRaw(S, K_abs, F);

  GN, mpN := RayClassGroup(frakN, [1 .. Degree(F)]);
  by_chi := AssociativeArray();
  for psi0 in S do
    psi := AssociatedPrimitiveCharacter(psi0);
    D := ZF !! (Norm(ZK!!Conductor(psi)) * rel_disc);
    if not IsIntegral(frakN/D) then
      continue;
    end if;
    B_N := <<mpN(g), psi(Integers(K_abs)!!mpN(g))^-1 * QuadraticCharacter(mpN(g), K)> : g in Generators(GN)>;
    chi_N := HeckeCharacter(frakN, [1 .. Degree(F)], B_N);
    if not IsDefined(by_chi, chi_N) then
      by_chi[chi_N] := [* *];
    end if;
    Append(~by_chi[chi_N], psi);
  end for;

  results := [* *];
  for chi_N -> psilist in by_chi do
    Mk_N := HMFSpace(GRing, frakN, [1,1], chi_N);
    new_forms := [* ThetaSeries(Mk_N, psi) : psi in psilist *];

    old_psis := PossibleGrossencharsOfRelQuadExt(K, frakN, [1,1], chi_N : AllowImprimitive:=true, GRing:=GRing);
    old_forms := [* *];
    for oldpsi in old_psis do
      Primitivize(oldpsi);
      Append(~old_forms, ThetaSeries(Mk_N, oldpsi));
    end for;

    matches := [];
    for nf in new_forms do
      found := false;
      for of in old_forms do
        if nf eq of then
          found := true;
          break;
        end if;
      end for;
      Append(~matches, found);
    end for;

    Append(~results, <chi_N, #new_forms, #old_forms, #new_forms eq #old_forms, matches, new_forms, old_forms>);
  end for;
  return results;
end intrinsic;

intrinsic RawDihedralCandidatesDebug(GRing::ModFrmHilDGRng, frakN::RngOrdIdl, Ntgt::RngIntElt) -> List
  {
    Debug/testing utility: like RawDihedralCandidates, but broken down per
    auxiliary K, and applies Basis() to the raw (pre-Inclusion) theta series
    produced for each K to check for linear redundancy among survivors that
    PrunedGrossencharsSetRaw (conjugate-pair dedup only) would not catch.
    Returns a List of <K, D_norms, num_raw_forms, num_after_basis>.
  }
  F := BaseField(GRing);
  ZF := Integers(F);
  Ks := QuadraticExtensionsWithConductor(frakN, [1 .. Degree(F)]);

  out := [* *];
  for K in Ks do
    S := TrivialDiagonalRestrictionGrossenchars(K, frakN, Ntgt);
    if #S eq 0 then continue; end if;
    K_abs := AbsoluteField(K);
    S := PrunedGrossencharsSetRaw(S, K_abs, F);

    ZK := Integers(K);
    rel_disc := Discriminant(ZK);
    // group surviving psi by D, tracking each psi's OWN independently
    // reconstructed chi_D (not shared/assumed equal across the group)
    psis_by_D := AssociativeArray();
    chiD_disagreements := 0;
    for psi0 in S do
      psi := AssociatedPrimitiveCharacter(psi0);
      D := ZF !! (Norm(ZK!!Conductor(psi)) * rel_disc);
      if not IsIntegral(frakN/D) then continue; end if;
      GD, mpD := RayClassGroup(D, [1 .. Degree(F)]);
      B_D := <<mpD(g), psi(Integers(K_abs)!!mpD(g))^-1 * QuadraticCharacter(mpD(g), K)> : g in Generators(GD)>;
      chi_D := HeckeCharacter(D, [1 .. Degree(F)], B_D);
      if not IsDefined(psis_by_D, D) then
        psis_by_D[D] := [* <psi, chi_D> *];
      else
        if chi_D ne psis_by_D[D][1][2] then
          chiD_disagreements +:= 1;
        end if;
        Append(~psis_by_D[D], <psi, chi_D>);
      end if;
    end for;

    num_raw := 0;
    num_after_basis := 0;
    D_norms := [];
    for D -> pclist in psis_by_D do
      num_raw +:= #pclist;
      Append(~D_norms, Norm(D));
      // build Mk_D/theta series per psi with ITS OWN chi_D, then group by
      // (as opposed to shared) chi_D before calling Basis
      by_chiD := AssociativeArray();
      for pc in pclist do
        psi, chi_D := Explode(pc);
        Mk_D := HMFSpace(GRing, D, [1,1], chi_D);
        f := ThetaSeries(Mk_D, psi);
        if not IsDefined(by_chiD, chi_D) then
          by_chiD[chi_D] := [];
        end if;
        Append(~by_chiD[chi_D], f);
      end for;
      for chi_D -> flist in by_chiD do
        if #flist gt 1 then
          num_after_basis +:= #Basis(flist);
        else
          num_after_basis +:= 1;
        end if;
      end for;
    end for;

    if num_raw eq 0 then continue; end if;
    if chiD_disagreements gt 0 then
      printf "  [K disc %o] WARNING: %o psi pairs at same D induced DIFFERENT chi_D\n",
        Discriminant(ZK), chiD_disagreements;
    end if;
    Append(~out, <K, D_norms, num_raw, num_after_basis>);
  end for;
  return out;
end intrinsic;
