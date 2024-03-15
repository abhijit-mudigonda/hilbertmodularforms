procedure test(K, L)
  // K - FldNum
  // L - FldNum

  K_embs := EmbeddingsCC(K);
  L_embs := EmbeddingsCC(L);

  for v in K_embs do
    extends_ct := 0;
    for w in L_embs do
      if Extends(w, v) then
        extends_ct +:= 1;
      end if;
    end for;
  end for;

  for F in [K, L] do
    assert #EmbeddingsCC(F) eq Degree(F);
    assert #[v : v in EmbeddingsCC(F) | IsReal(v)] eq #RealPlaces(F);

    for f in EmbeddingsCC(F) do
      assert (Conjugate(f)`Place eq f`Place);
      assert #[g : g in EmbeddingsCC(F) | IsConjugate(f, g)] eq 1;
      a := PrimitiveElement(F);
      if IsReal(f) then
        assert f(a) eq Evaluate(a, f`Place);
      else
        assert (Conjugate(f)`Conj ne f`Conj);
        b := Evaluate(a, f`Place);
        assert f(a) eq ((f`Conj) select Conjugate(b) else b);
        assert Conjugate(f)(a) eq Conjugate(f(a));
      end if;
    end for;
  end for;
end procedure;

K := QuadraticField(5);
L := CyclotomicField(5);
test(K, L);

R<x> := PolynomialRing(Rationals());
K := NumberField(x^3-2);
S<y> := PolynomialRing(K);
L_rel := ext<K | y^3 - (3+K.1-K.1^2)>;
L := AbsoluteField(L_rel);
assert not IsGalois(L);
test(K, L);
