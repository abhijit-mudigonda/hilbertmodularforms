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
