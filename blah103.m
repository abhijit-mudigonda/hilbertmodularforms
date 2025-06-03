load "config.m";

import "ModFrmHil/level.m" : InducedH1Internal, CompleteRelationFromUnit;
import "ModFrmHil/weight_rep.m" : weight_rep_dim;

SetVerbose("ModFrmHil", 3);

// Given ideal data X and X_m, adjusts the coset reps of X
// so that the coset reps of X_m are a subset of those of X.
function update_coset_reps(X, X_m)
  assert X`Ideal eq Discriminant(X_m`QuaternionOrder) * X_m`Ideal;
  O := X`QuaternionOrder;
  O_m := X_m`QuaternionOrder;
  ZF := Integers(BaseField(Algebra(O)));
  subgp_coset_idxs := [];
  for w in X_m`P1Elements do
    // print "hoiiiiiiiiiiiiiiiiiiiiiiiiiiiiii";
    // print "w", w;
    _, c := Explode(CosetRepsByP1(X_m)[w]);
    // _, wp := X_m`P1Rep(X_m`ResidueMap(c)[2], false, false);
    // assert wp eq w;
    // assert Parent(c) eq O_m;
    // print "X_m`ResidueMap(c)[2]", X_m`ResidueMap(c)[2];
    // print "X`ResidueMap(c)[2]", X`ResidueMap(c)[2];
    // print "c", c;
    // assert &+[Eltseq(c)[i] * O_m.i : i in [1 .. 4]] eq c;

    c := O!c;
    // print "O!c", c;
    u := X`ResidueMap(c)[2];
    // print "u", u;
    assert ideal<ZF | u[1], u[2]> eq 1*ZF;
    _, v := X`P1Rep(u, false, false);
    // print "v", v;
    assert v in Keys(X`CosetRepsByP1);
    j := Index(X`P1Elements, v);
    // assert not (j in subgp_coset_idxs);
    assert X`CosetReps[j] eq X`CosetRepsByP1[v][2];
    // set the coset rep at v to c 
    X`CosetRepsByP1[v] := <j, c>;
    X`CosetReps[j] := c;
    Append(~subgp_coset_idxs, j);
  end for;

  print subgp_coset_idxs;

  for c in X_m`CosetReps do
    idx := Index(X`CosetReps, O!c);
    assert idx ne 0;
    assert idx in subgp_coset_idxs;
  end for;

  return subgp_coset_idxs;
end function;

function shapiro_matrix(X, X_m, k, subgp_coset_idxs)
  Htilde, mH := InducedH1Internal(X, k);
  Htilde_m, mH_m := InducedH1Internal(X_m, k);
  print "Htilde", Nrows(Htilde), Ncols(Htilde);
  print "Htilde_M", Nrows(Htilde_m), Ncols(Htilde_m);

  Z := Domain(mH);
  Z_m := Domain(mH_m);

  U, _, m := Group(X`FuchsianGroup);
  U_m, _, m_m := Group(X_m`FuchsianGroup);

  dim_W := weight_rep_dim(k);
  dim_indW_X := dim_W * #X`P1Elements;
  dim_indW_Xm := dim_W * #X_m`P1Elements;

  // The Shapiro isomorphism is a restriction (to the smaller group)
  // followed by a projection (from the induced module to the original module)
  O_m := X_m`QuaternionOrder;
  O := X`QuaternionOrder;
  assert &and[m_m(U_m.i) in O_m : i in [1 .. #Generators(U_m)]];
  assert &and[IsCoercible(O, m_m(U_m.i)) : i in [1 .. #Generators(U_m)]];
  rel_mtrxs := [CompleteRelationFromUnit(X, m_m(U_m.i), k) : i in [1 .. #Generators(U_m)]];
  res_mtrx := HorizontalJoin(rel_mtrxs);
  assert Nrows(res_mtrx) eq Ncols(res_mtrx);
  assert Nrows(res_mtrx) eq #Generators(U_m) * dim_indW_X;

  print "dim_indW_X", dim_indW_X;
  subgp_mtrx_idxs := &cat[&cat[[(j * dim_indW_X + ((i - 1) * dim_W) + 1) .. (j * dim_indW_X + (i * dim_W))] : i in subgp_coset_idxs] : j in [0, 1]];
  assert #SequenceToSet(subgp_mtrx_idxs) eq #subgp_mtrx_idxs;
  print "subgp_mtrx_idxs", subgp_mtrx_idxs;

  res_mtrx := Submatrix(res_mtrx, [1 .. Nrows(res_mtrx)], subgp_mtrx_idxs);
  assert Ncols(res_mtrx) eq #Generators(U_m) * dim_indW_Xm;

  // The rows of Htilde are lifts of a basis of Z/B to Z.
  // res_mtrx takes Z to Z_m
  A := Htilde * res_mtrx;
  S := Matrix([mH_m(A[i]) : i in [1 .. Nrows(A)]]);

  return S;
end function;

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

B := QuaternionAlgebra<F | -3, -12*F.1>;
O := MaximalOrder(B);
Gamma := FuchsianGroup(O);
U, mm, m := Group(Gamma);
N := 6*ZF;
k := [2,2,4];
chi := HeckeCharacterGroup(N, [1,2,3]).0;
X := cIdealDatum(Gamma, N);

O_3 := Order(O, 3*ZF);
Gamma_3 := FuchsianGroup(O_3);
U_m, mm_m, m_m := Group(Gamma_3);
H_3 := HeckeCharacterGroup(2*ZF, [1,2,3]);
chi_3 := Restrict(chi, H_3);
X_3 := cIdealDatum(Gamma_3, 2*ZF : residue_map:=X`ResidueMap);

subgp_coset_idxs := update_coset_reps(X, X_3);

pp := Factorization(7*ZF)[1][1];
T7 := HeckeMatrix2(Gamma_3, 2*ZF, pp, k, chi_3);

pp := Factorization(13*ZF)[1][1];
T13 := HeckeMatrix2(Gamma_3, 2*ZF, pp, k, chi_3);

Tinf := HeckeMatrix2(Gamma_3, 2*ZF, "Infinity", k, chi_3);

pp := 2*ZF;
// T2 := HeckeMatrix2(Gamma_3, 2*ZF, pp, k, chi_3);

pp := 3*ZF;
T3 := HeckeMatrix2(Gamma, 6*ZF, pp, k, chi);

S := shapiro_matrix(X, X_3, k, subgp_coset_idxs);
