declare attributes AlgQuat : Splittings;

function GetOrMakeP1_new(Gamma, N)
  // Gamma - GrpPSL2
  // N - RngOrdIdl
  //
  // Returns the cached output of ProjectiveLine(Gamma, N)
  Z_F := Order(N);
  Z_FN := quo<Z_F | N>;
  if not assigned Gamma`P1s_new then
    Gamma`P1s_new := AssociativeArray();
  end if;
  if IsDefined(Gamma`P1s_new, N) then
    return Explode(Gamma`P1s_new[N]);
  else
    P1N, P1Nrep := ProjectiveLine(Z_FN);
    Gamma`P1s_new[N] := <P1N, P1Nrep>;
    return P1N, P1Nrep;
  end if;
end function;

function FastSymmetricPower(mat, n)
  case n:
    when 0:
      return MatrixRing(BaseRing(mat), 1)!1;
    when 1:
      return mat;
    else
      return SymmetricPower2(mat, n);
  end case;
end function;

function weight_map_arch(b, m, n)
  // b - AlgQuatElt or AlgAssVOrdElt
  // m - SeqEnum[FldRatElt] or SeqEnum[RngIntElt]
  // n - SeqEnum[RngIntElt]
  //
  // returns a matrix corresponding to the action of b on
  // Wk = ⊗_i (Sym^(n_i) ⊗ det^(m_i)),
  // where b acts on the ith component via its image under
  // the ith embedding B -> M_2(K) coming from Splittings(B).
  //
  // In practice, n = k - 2 and m depends on the choice of central character.
  // m will usually be integral but in the non-paritious case it is not possible
  // for it to be integral (but even then it will be at worst half-integral). 
  //
  // Some of the calls to this function are for use as the coefficient module 
  // for group cohomology of an arithmetic Fuchsian group.
  // In this setting, the reduced norm is always 1 so we can ignore the determinant terms. 
  //
  // TODO abhijitm - the nebentypus will also get added here. 
  d := #m;
  splittings, K, _ := Splittings(Parent(b));
  M := MatrixRing(K, 1)!1;
  for l := d to 1 by -1 do 
    // this casework looks gross because this function needs to be performant. 
    if m[l] eq 0 and n[l] eq 0 then
      // don't need to modify M
      continue;
    else
      mat := splittings[l](b);
      if n[l] eq 0 then
        Ml := MatrixRing(K, 1)!(Determinant(mat)^m[l]);
      elif m[l] eq 0 then
        Ml := FastSymmetricPower(mat, n[l]);
      else
        // TODO abhijitm m[l] is a fraction then this probably needs to be in a
        // bigger number field, but I'm going to ignore this issue for now. 
        Ml := Determinant(mat)^m[l] * FastSymmetricPower(mat, n[l]);
      end if;
      // Tensor product is associative; for efficiency always do
      // TensorProduct(small_mat, large_mat)
      M := TensorProduct(Ml, M);
    end if;
  end for;
  return M;
end function;

//-------------
//
// Right action functions.
//
//-------------

// Gamma`LevelRPAs_new stores a dictionary keyed by level N.
// Gamma`LevelRPAs_new[N] contains a dictionary keyed by 
// the integers [-#gens ..., -1, 1, ..., #gens], where gens is
// Generators(Group(Gamma)).
//
// At a positive integer i, the dictionary stores the permutation matrix 
// induced by right action of the ith generator on coset representatives
// (which should be something like Gamma(N) \ Gamma). 
function RightPermutationActions(X)
  // X::IdealDatum
  Gamma := X`FuchsianGroup;
  N := X`Ideal;
  if not assigned Gamma`LevelRPAs_new then
    Gamma`LevelRPAs_new := AssociativeArray();
  end if;

  if IsDefined(Gamma`LevelRPAs_new, N) then
    return Gamma`LevelRPAs_new[N];
  end if;

  vprintf ModFrmHil: "Computing right permutation actions .................. ";
  time0 := Cputime();

  U, m := Group(Gamma);
  RPAs := AssociativeArray();
  for i := 1 to #Generators(U) do
    delta := Quaternion(m(U.i));
    perm := [];
    for alphai in X`CosetReps do
      _, v := X`P1Rep(X`ResidueMap(alphai*delta)[2], false, false);
      Append(~perm, Index(X`P1Elements, v));
    end for;
    RPAs[i] := PermutationSparseMatrix(Integers(), SymmetricGroup(#X`P1Elements)!perm);
    RPAs[-i] := PermutationSparseMatrix(Integers(), SymmetricGroup(#X`P1Elements)!perm^-1);
  end for;

  vprintf ModFrmHil: "Time: %o\n", Cputime(time0);

  Gamma`LevelRPAs_new[N] := RPAs;
  return RPAs;
end function;

intrinsic Splittings(O::AlgAssVOrd) -> SeqEnum[Map], FldNum, FldNum
  {}
  return Splittings(Algebra(O));
end intrinsic;

intrinsic Splittings(B::AlgQuat) -> SeqEnum[Map], FldNum, FldNum
  {
    input: 
      B - A quaternion algebra defined over a degree n field F.
    returns:
      A SeqEnum of n maps from B to M2(K), where K is the minimal field
      over which such a sequence can be defined. Each map in the sequence
      corresponds to one of the infinite places of F. 

      We also return the field K as well as (TODO abhijitm no clue why but 
      I didn't want to change it) an unoptimized version of K. 
  }
  if assigned B`Splittings then
    return Explode(B`Splittings);
  end if;

  F := BaseField(B);
  // define weight_base_field = extension K/F containing Galois closure of F and 
  // containing a root of every conjugate of the minimal polynomial of B.1
  if assigned F`SplittingField then
    K,rts := Explode(F`SplittingField);
  else
    K,rts := SplittingField(F : Abs := true, Opt := false);
    F`SplittingField := <K, rts>;
  end if;
  embeddings_F_to_K := [hom<F->K | r> : r in rts];
  B1coeffs := Coefficients(MinimalPolynomial(B.1));
  alphas := [K| ];
  for FtoK in embeddings_F_to_K do
    hh := PolynomialRing(K)! [c@FtoK : c in B1coeffs];
    if IsIrreducible(hh) then
       K := ext<K|hh>;
       alphas := ChangeUniverse(alphas,K) cat [K.1];
    else
       Append(~alphas, Roots(hh)[1][1]);
    end if;
  end for;
  // make weight_base_field an (optimized) absolute field, for efficiency in later calculations 
  weight_field := K; // names appears in verbose output 
                     // TODO abhijitm I don't really get why we can't let K be the optimized
                     // representation but I'm too lazy to think about it now and this is how
                     // it was so we return this separately. 
  K := AbsoluteField(K);
  K := OptimizedRepresentation(K);
  embeddings_F_to_K  :=  [hom<F->K | K!r> : r in rts]; // same embeddings, now into extended field K

  require B.1*B.2 eq B.3 : "We assume B.1 * B.2 == B.3 when defining\
    the splitting homomorphisms.";
  splitting_seq := [];
  for i := 1 to Degree(F) do
    h := embeddings_F_to_K[i];
    // need a splitting homomorphism (B tensor K) -> Mat_2(K) whose restriction to K is h 
    alpha := alphas[i];
    b := K! h(F!(B.2^2));
    iK := Matrix(K, 2, [alpha, 0, 0, -alpha]); 
    jK := Matrix(K, 2, [0, 1, b, 0]); 
    kK := iK*jK;
    assert K! h(B.3^2) eq (kK^2)[1,1]; 
    Append(~splitting_seq, 
        map< B -> MatrixRing(K,2)|
          q:-> h(s[1])+h(s[2])*iK+h(s[3])*jK+h(s[4])*kK where s := Eltseq(q) >);
  end for;
  B`Splittings := <splitting_seq, K, weight_field>;
  return splitting_seq, K, weight_field;
end intrinsic;
