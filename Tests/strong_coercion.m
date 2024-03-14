CC_THRESHOLD := 10^-10;

procedure test(a, K1, K2 : v1:=DefaultMarkedEmbedding(K1), v2:=DefaultMarkedEmbedding(K2))
  // a::FldElt - Element of K1
  // K1::Fld  
  // K2::Fld 
  // v1::EmbedNumElt
  // v2::EmbedNumElt
  assert a in K1;
  able, b := IsStrongCoercible(K2, a : v:=v1, w:=v2); 
  if able then
    a_eval := (K1 ne Rationals()) select v1(a) else a;
    b_eval := (K2 ne Rationals()) select v2(b) else b;
    assert Abs(a_eval - b_eval) lt CC_THRESHOLD;
    c := StrongCoerce(K1, b : v:=v2, w:=v1);
    assert c eq a;
  end if;
end procedure;

Q := RationalField();
F := QuadraticField(5);
R<x> := PolynomialRing(Rationals());
K := NumberField(x^3 - x^2 - 2*x + 1);
L := CyclotomicField(5);
H := Compositum(F, L);
M := Compositum(K, L);

// FldRat <-> FldRat
test(163/1729, Q, Q);

// FldRat <-> FldQuad
test(163/1729, Q, F);

// FldRat <-> FldCyc
test(163/1729, Q, L);

// FldQuad <-> FldNum
test(163 + 1729*F.1, F, H);

// FldQuad <-> FldCyc
// TODO abhijitm we omit this test for now because if x is in K
// and L is a cyclotomic field containing K, then
// L!x will succeed but K!(L!x) will fail

// FldCyc <-> FldNum
test(16 + 3*L.1 + 17*L.1^2 + 2*L.1^3 + 9*L.1^4, L, H);

// FldNum <-> FldNum
test(163 + 17*K.1 + 29*K.1^2, K, M);

// StrongMultiply
assert StrongMultiply([* K.1, L.1^3, K.1^-1, L.1^-3, 18/41 *] : K:=M) eq 18/41;

v_K := DefaultMarkedEmbedding(K);
v_L := DefaultMarkedEmbedding(L);
w := DefaultMarkedEmbedding(M);
a := K.1;
b := L.1;
c := StrongMultiply([* a, b *] : K:=M);
assert Abs((c @ w) - (a @ v_K) * (b @ v_L)) lt CC_THRESHOLD;

B := ListToStrongCoercedSeq([* 1, 2/3, K.1, L.1 *]);
assert IsIsomorphic(Parent(B[1]), Compositum(K, L));

// non-Galois <-> Galois
R<x> := PolynomialRing(Rationals());
K := NumberField(x^3-2);
L := SplittingField(K);
vs := EmbeddingsCC(K);
w := DefaultMarkedEmbedding(L);
for v in vs do
  test(K.1 - 3, K, L : v1:=v, v2:=w);
end for;

// Galois <-> non-Galois
F := QuadraticField(5);
R<x> := PolynomialRing(F);
K_rel := ext<F | x^2 + 1/22*(3*F.1 + 23)>;
K := AbsoluteField(K_rel);
test(F.1, F, K : v2:=EmbeddingsCC(K)[2]);

// non-Galois <-> non-Galois
R<x> := PolynomialRing(Rationals());
K := NumberField(x^3-2);
S<y> := PolynomialRing(K);
L_rel := ext<K | y^3 - (3+K.1-K.1^2)>;
L := AbsoluteField(L_rel);
assert not IsGalois(L);

for v in EmbeddingsCC(K) do
  for w in EmbeddingsCC(L) do
    test(K.1, K, L : v1:=v, v2:=w);
    assert IsStrongCoercible(L, K.1 : v:=v, w:=w) eq Extends(w, v); 
  end for;
end for;
