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
  return splitting_seq, K, weight_field;
end intrinsic;
