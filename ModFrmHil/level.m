freeze;

//////////////////////////////////////////////////////////////////////////////
//
// Hecke operators with level structure for Indefinite Algorithm
// January 2009, John Voight
//
//////////////////////////////////////////////////////////////////////////////

import !"Geometry/ModFrmHil/proj1.m" : residue_class_reps;
import !"Geometry/ModFrmHil/indefinite.m" : ElementOfNormMinusOne, LeftIdealGens;
import !"Geometry/ModFrmHil/hecke.m" : pseudo_inverse;

import !"Geometry/ModFrmHil/level.m" : FindGammas;

import !"Algebra/AlgQuat/enumerate.m" :
             EnumerativeSearchInternal, ReducedBasisInternal;
import !"Geometry/GrpPSL2/GrpPSL2Shim/domain.m" : Vertices;
import "weight_rep.m" : GetOrMakeP1_new, Gamma0Cosets, RightPermutationActions;

ConjugationPermutationActions := function(Gamma, N, Z_FN, iota, P1N, cosets, P1Nrep);
  if not assigned Gamma`LevelCPAs_new then
    Gamma`LevelCPAs_new := AssociativeArray();
  end if;

  if IsDefined(Gamma`LevelCPAs_new, N) then
    return Explode(Gamma`LevelCPAs_new[N]);
  end if;

  Z_F := BaseRing(BaseRing(Gamma));
  bas, n_seq := residue_class_reps(N);
  Rset:=[[s[m]: m in [1..#s]]: s in Set(CartesianProduct(<[0..n_seq[l]-1]: l in [1..#n_seq]>))];

  iotaalphavs := [];
  for alphai in cosets do
    _, v := P1Nrep(iota(alphai)[2], false, false);
    Append(~iotaalphavs, [Z_F!t : t in Eltseq(v)]);
  end for;

  qcnt := 0;
  CPAs1bas := [];
  for q in bas do
    qcnt +:= 1;
    perm := [];
    for w in iotaalphavs do
      _, v := P1Nrep([w[1], w[1]*q + w[2]], false, false);
      Append(~perm, Index(P1N, v));
    end for;
    perm := SymmetricGroup(#cosets)!perm;
    Append(~CPAs1bas, perm);
  end for;

  Z_FNstar, mZ_FNstar := UnitGroup(Z_FN);
  basmult := [Z_F!mZ_FNstar(Z_FNstar.i) : i in [1..#Generators(Z_FNstar)]];
  qcnt := 0;
  CPAs2bas := [];
  for q in basmult do
    qcnt +:= 1;
    perm := [];
    for w in iotaalphavs do
      _, v := P1Nrep([w[1], q*w[2]], false, false);
      Append(~perm, Index(P1N, v));
    end for;
    perm := SymmetricGroup(#cosets)!perm;
    Append(~CPAs2bas, perm);
  end for;
  
  Q1 := [Z_FN!x : x in Rset];
  CPAs1 := [];
  for i := 1 to #Rset do
    perm := &*[CPAs1bas[j]^Rset[i][j] : j in [1..#CPAs1bas]];
    Append(~CPAs1, perm);
  end for;
  ChangeUniverse(~Q1, Z_FN);

  Q2 := [];
  CPAs2 := [];
  for i := 1 to #Rset do
    z := Z_FN!Rset[i];
    if IsUnit(z) then
      perm := &*[CPAs2bas[j]^zseq[j] : j in [1..#CPAs2bas]] where zseq is Eltseq(z@@mZ_FNstar);
      Append(~CPAs2, perm);
      Append(~Q2, z);
    end if;
  end for;

  Gamma`LevelCPAs_new[N] := <Q1, CPAs1, Q2, CPAs2>;
  return Q1, CPAs1, Q2, CPAs2;  
end function;

//-------------
//
// Interface with fundamental domain algorithm.
//
//-------------

InducedRelation := function(rel, RPAs, RPAsinv : IsTrivialCoefficientModule:=false);
  rpa1 := RPAs[1];
  R    := BaseRing(rpa1);
  nr   := Nrows(rpa1);
  mats := [SparseMatrix(R, nr, nr) : i in [1..#RPAs]];

  I    := IdentitySparseMatrix(R, nr);
  if IsTrivialCoefficientModule then
    for i := 1 to #rel do
      absi := Abs(rel[i]);
      if rel[i] lt 0 then
        mats[absi] -:= I;
      else
        mats[absi] +:= I;
      end if;
    end for;
  else
    g := I;
    for i := #rel to 1 by -1 do
      absi := Abs(rel[i]);
      if rel[i] lt 0 then
        mats[absi] -:= RPAsinv[absi]*g;
      else
        mats[absi] +:= g;
      end if;
      if Sign(rel[i]) eq 1 then
        g := RPAs[absi]*g;
      else
        g := RPAsinv[absi]*g;
      end if;
    end for;
  end if;
  return VerticalJoin(mats), rel;
end function;

CompleteRelationFromUnit := function(Gamma, alpha, RPAs, RPAsinv : IsTrivialCoefficientModule:=false);
  reldata := ShimuraReduceUnit(Gamma!alpha);
  assert IsScalar(Quaternion(reldata[1]));
  rel := reldata[3];

  mat, rel := InducedRelation(rel, RPAs, RPAsinv : IsTrivialCoefficientModule:=IsTrivialCoefficientModule);
  return mat, rel;
end function;

//-------------
//
// Cohomology module.
//
//-------------

function InducedH1Internal(Gamma, N, cosets, RPAs, RPAsinv);
  if assigned Gamma`LevelH1s then
    for H1 in Gamma`LevelH1s do
      if H1[1] eq N then
        return H1[2], H1[3];
      end if;
    end for;
  end if;

  U, m := Group(Gamma);
  d := #Generators(U);
  gammagens := [Quaternion(m(U.i)) : i in [1..d]];

  R := HorizontalJoin(
    [InducedRelation(Eltseq(LHS(rel)), RPAs, RPAsinv) : rel in Relations(U)]);
  Z := Kernel(R);

  I := IdentitySparseMatrix(BaseRing(RPAs[1]), Nrows(RPAs[1]));
  coB := HorizontalJoin([r - I : r in RPAs]);
  coB := [Z!coB[i] : i in [1..Nrows(coB)]];
  ZcoB := sub<Z | coB>;

  H, mH := quo<Z | ZcoB>;

  ZB := Basis(Z);
  S := EchelonForm(Matrix(coB));
  RemoveZeroRows(~S);
  piv := [Min(Support(S[i])): i in [1..Nrows(S)]];
  assert #SequenceToSet(piv) eq #piv;
  Htilde := [ZB[i] : i in [1..#ZB] | 
                        Min(Support(ZB[i])) notin piv];
  if #Htilde gt 0 then
    mHtilde := Matrix([mH(h) : h in Htilde]);
    assert Abs(Determinant(mHtilde)) eq 1;
    Htilde := mHtilde^(-1)*Matrix(Htilde);
    Htilde := [Htilde[i] : i in [1..Nrows(Htilde)]];
  end if;

  if assigned Gamma`LevelH1s then
    Append(~Gamma`LevelH1s, <N, Htilde, mH>);
  else
    Gamma`LevelH1s := <<N, Htilde, mH>>;
  end if;

  return Htilde, mH;
end function;

function InducedH1(Gamma, Gammap, N, cosets, cosetsp, RPAs, RPAsinv, RPAsp, RPAspinv);
  Htilde, mH := InducedH1Internal(Gamma, N, cosets, RPAs, RPAsinv);
  Htildep, mHp := InducedH1Internal(Gammap, N, cosetsp, RPAsp, RPAspinv);

  return Htildep, mH;
end function;

//-------------
//
// Main loop.
//
//-------------

HeckeMatrix1 := function(O_mother, N, ell, ind, indp, ridsbasis, iotaell : ellAL := false);
  // Initialization.
  Gamma_mother := O_mother`FuchsianGroup;
  assert O_mother`RightIdealClasses[ridsbasis][4];
  rids := O_mother`RightIdealClasses[ridsbasis];

// GetMemoryUsage(); MemProfile();

  P1N, P1Nrep := GetOrMakeP1_new(Gamma_mother, N);

  B := Algebra(O_mother);
  F := BaseRing(B);
  Z_F := MaximalOrder(F);
  Z_FN := quo<Z_F | N>;

  O := rids[3][ind];
  Op := rids[3][indp];
  Gamma := O`FuchsianGroup;
  Gammap := Op`FuchsianGroup;

  J := rids[2][ind];
  Jp := rids[2][indp];

  FeqQQ := F cmpeq Rationals();

  if FeqQQ then 
    JJp := Jp;
  else
    JJp := Jp*LeftInverse(J);
  end if;

  elleqoo := ell cmpeq "Infinity";
  ellU := not elleqoo and ell + Discriminant(O)/Discriminant(B)*N eq ell;
  inNormSupport := not elleqoo and (iotaell cmpeq [] or 
    Gcd(Integers()!Norm(ell),Integers()!Norm(rids[1][ind]*rids[1][indp])) ne 1);

  U, _, m := Group(Gamma);

  Uside := Gamma`ShimGroupSidepairs;
  mside := Gamma`ShimGroupSidepairsMap;
  n := #Generators(U);
  lifts := [m(U.i) : i in [1..n]];

  IsLevelOne := Norm(N) eq 1;

  // Check or precompute level structure.
  _, iota := ResidueMatrixRing(O, N);
  _, iotap := ResidueMatrixRing(Op, N);
  Z_FN := quo<Z_F | N>;

  cosets := Gamma0Cosets(Gamma, N, Z_FN, iota, P1N, P1Nrep);
  cosetsp := Gamma0Cosets(Gammap, N, Z_FN, iotap, P1N, P1Nrep);

  RPAs, RPAsinv := RightPermutationActions(Gamma, N, Z_FN, iota, P1N, cosets, P1Nrep);
  RPAsp, RPAspinv := RightPermutationActions(Gammap, N, Z_FN, iotap, P1N, cosetsp, P1Nrep);

  D := Parent(Gamma`ShimFDDisc[1]); 

  // There are two methods to compute the Hecke operator.
  // One works in the situation when ell is coprime to N and the support of the
  // right ideal classes; it runs zippily.
  // In the other situations, one must "work hard", which means unpacking 
  // Shapiro's lemma and dealing with many issues of normalizations.

  // If ell = oo, we need to decide if we have to work hard or not.
  if elleqoo then
    P1ell := [Infinity()];
    numP1 := 1;
    ellcosets := [B!1];

    if FeqQQ then
      alphap := ElementOfNormMinusOne(O);
    else
      _, alphap := IsPrincipal(JJp, Gammap : Strict := false);
    end if;

    if not FeqQQ and not -1 in RealSigns(Norm(alphap)) then
      assert ind eq indp;
      alphap := ElementOfNormMinusOne(O);
      assert alphap in O and IsUnit(Z_F!Norm(alphap)) and
             -1 in RealSigns(Norm(alphap));
    end if;

    // Ensure alpha is trivial at N.
    iotaalphap := iotap(alphap);
    v := [iotaalphap[2,1], -iotaalphap[1,1]];
    if not IsLevelOne then
      _, v := P1Nrep(v, false, false);
      c := cosetsp[Index(P1N, v)];
      alphap := c*alphap;
    end if;

    lambda := alphap;
    lambdas := [lambda];
    NNlambda := Norm(Norm(lambda));
    NNlambda := Numerator(NNlambda)*Denominator(NNlambda);

    ooinNormSupport := Gcd(NNlambda,Integers()!Norm(rids[1][ind]*rids[1][indp])) ne 1;
  else
    ooinNormSupport := false;
  end if;

  // Catch the cases where we work hard:
  // (1) ell = oo and the representative element of negative norm is not coprime;
  // (2) ell divides N (but not ell divides D)--this includes the Hecke operator case
  //     ellU and the Atkin-Lehner case ellAL.
  if (elleqoo and ooinNormSupport) or (ellAL and Valuation(N,ell) gt 0) or ellU then
    if elleqoo then
      numP1 := 1;
    elif ellAL then
      numP1 := 1;
    elif ellU then
      numP1 := Norm(ell);
    end if;

    // Work from the definition to find lambdas.
    // If ell = oo, we already have lambda.  
    if ellAL then
      ee := Valuation(Discriminant(O)*N, ell);
      lambda := LeftIdealGens(Gamma, ell, JJp, 1, O, Op, iotaell : Slow := true, ALval := ee)[1];

      lambdas := [];
      _, iotapNell := ResidueMatrixRing(O, N/ell^ee);
      for c := 1 to #cosetsp do
        if Valuation(iotap(cosetsp[c]*lambda)[2,1],ell) ge ee and
           Valuation(iotap(cosetsp[c]*lambda)[2,2],ell) ge ee and
           iotapNell(cosetsp[c]*lambda)[2,1] in N/ell^ee then
          Append(~lambdas, cosetsp[c]*lambda);
        end if;
      end for;
      assert #lambdas eq 1;
    elif ellU then
      Z_Fell := quo<Z_F | ell>;
      P1ellfull, P1ellfullrep := GetOrMakeP1_new(Gamma_mother, ell);
      ellcosetsfull := Gamma0Cosets(Gamma, ell, Z_Fell, iotaell, P1ellfull, P1ellfullrep);
      ooind := 1;
      while ooind le #P1ellfull do
        if Valuation(P1ellfull[ooind][1],ell) gt 0 then // Should be infinity, in fact!
          break;
        end if;
      end while;
      P1ell := [P1ellfull[c] : c in [1..#P1ellfull] | c ne ooind];
      ellcosets := [ellcosetsfull[c] : c in [1..#P1ellfull] | c ne ooind];

      lambda := LeftIdealGens(Gamma, ell, JJp, 1, O, Op, iotaell);
      lambdas := [lambda*ellcosets[c] : c in [1..numP1]];

      for i := 1 to #lambdas do
        lambdap := lambdas[i];
        iotalambdap := iotap(lambdap);
        v := [iotalambdap[2,1], -iotalambdap[1,1]];
        _, v := P1Nrep(v, false, false);
        c := cosetsp[Index(P1N, v)];
        lambdas[i] := c*lambdas[i];
      end for;
    end if;

    Y_U := [];

// GetMemoryUsage(); MemProfile();

    vprintf ModFrmHil: "Computing operator the hard way ...................... ";
    vtime ModFrmHil:

    for i in [1..n] do
      G := [];
      for k in [1..#cosets] do
        Gk := [];
        iotalik := iota(lifts[i]*cosets[k]^(-1));
        v := [iotalik[2,1], -iotalik[1,1]];
        _, v := P1Nrep(v, false, false);
        ci := Index(P1N, v);
        liftsik := O!(cosets[ci]*lifts[i]*cosets[k]^(-1));
        Y_Opi := [];

        for j in [1..numP1] do
          if elleqoo or ellAL then
            c := 1;
          else
            iotadelta := iotaell(lambdas[j]*liftsik);
            bl, v := P1ellfullrep(iotadelta[1], false, false);
            c := Index(P1ell, v);
          end if;
          y := Op!(lambdas[j]*liftsik*lambdas[c]^(-1));
          y,rel := CompleteRelationFromUnit(Gammap, y, RPAsp, RPAspinv : IsTrivialCoefficientModule := false);
          Append(~Gk, y);
        end for;

        Append(~G, Gk);
      end for;
      Append(~Y_U, G);
    end for;

    Htilde, mH := InducedH1(Gamma, Gammap, N, cosets, cosetsp, RPAs, RPAsinv, RPAsp, RPAspinv);

    if #Htilde eq 0 then
      return [];
    else
      // There is some normalization issue that I'm missing.  It should work just taking the 
      // first column submatrix, but fails in some isolated cases.  Shapiro's lemma
      // works acceptably for any choice of constant element, so once we find one that 
      // works, we should be OK...
      for t := 1 to Ncols(Y_U[1][1][1]) do 
        M := HorizontalJoin([ HorizontalJoin([ &+[ ColumnSubmatrix(Y_U[i][k][j],t,1) : j in [1..numP1]] : k in [1..#cosets]]) : i in [1..n] ]);
        MH := Matrix(Htilde)*M;
        if &and[MH[i] in Domain(mH) : i in [1..#Htilde]] then
          MM := Matrix([mH(MH[i]) : i in [1..#Htilde]]);
          return MM;
        end if;
      end for;
      error "No column submatrix worked!?  This is a serious error; please report.";
    end if;
  end if;

// GetMemoryUsage(); MemProfile();

  // We've ruled out some "work hard" cases; we'll try to use as much optimization as possible.
  // We still have to do some extra computing if ell is in the support of the ideal classes.
  // First, if ell <> oo, we need to get our lambdas.
  if not elleqoo then
    if ellAL then
      numP1 := 1;
    else
      numP1 := Norm(ell)+1;
    end if;

    if ellAL or inNormSupport then
      if ellAL then
        // We've covered the case when ell divides N, so this is only the case ell divides D;
        // hence ell exactly divides D.
        lambdas := LeftIdealGens(Gamma, ell, JJp, 1, O, Op, iotaell : Slow := true, ALval := 1);
      elif inNormSupport then
        lambdas := LeftIdealGens(Gamma, ell, JJp, 1, O, Op, iotaell : Slow := true);
      end if;

      // Ensure lambda is trivial at N.
      if not IsLevelOne then
        for i := 1 to #lambdas do
          lambdap := lambdas[i];
          iotalambdap := iotap(lambdap);
          v := [iotalambdap[2,1], -iotalambdap[1,1]];
          _, v := P1Nrep(v, false, false);
          c := cosetsp[Index(P1N, v)];
          lambdas[i] := c*lambdas[i];
        end for;
      end if;
    else // Go for the fast code.
      lambda := LeftIdealGens(Gamma, ell, JJp, 1, O, Op, iotaell);
      P1ell, P1ellrep := GetOrMakeP1_new(Gamma_mother, ell);
      Z_Fell := quo<Z_F | ell>;
      ellcosets := Gamma0Cosets(Gamma, ell, Z_Fell, iotaell, P1ell, P1ellrep);

      // Ensure lambda is trivial at N.
      if not IsLevelOne then
        levelmults := [];
        for i := 1 to #ellcosets do
          lambdap := lambda*ellcosets[i];
          iotalambdap := iotap(lambdap);
          v := [iotalambdap[2,1], -iotalambdap[1,1]];
          _, v := P1Nrep(v, false, false);
          c := cosetsp[Index(P1N, v)];
          Append(~levelmults, c);
        end for;
      else
        levelmults := [Op!1 : i in [1..#ellcosets]];
      end if;
    end if;
  else
    levelmults := [Op!1];
  end if;

// GetMemoryUsage(); MemProfile();

  vprintf ModFrmHil: "Computing conjugation actions ........................ ";
  vtime ModFrmHil:
  if not IsLevelOne then
    Q1, CPAs1, Q2, CPAs2 := ConjugationPermutationActions(Gammap, N, Z_FN, iotap, P1N, cosetsp, P1Nrep);

    Zp := [];
    for j in [1..numP1] do
      if inNormSupport then 
        iotaj := iotap(lambdas[j]);
      else
        iotaj := iotap(levelmults[j]*lambda*ellcosets[j]);
      end if;
      xinv := (Z_FN!iotaj[1,1])^(-1);
      perm1 := CPAs1[Index(Q1, Z_FN!iotaj[1,2]*xinv)];
      perm2 := CPAs2[Index(Q2, Z_FN!iotaj[2,2]*xinv)];
      Append(~Zp, PermutationSparseMatrix(Integers(), perm2 * perm1));
    end for;
  end if;

// GetMemoryUsage(); MemProfile();

  Y_Op := [];
  X := [];
  vprintf ModFrmHil: "Defining maps for relations from units ............... ";
  vtime ModFrmHil:
  for i in [1..n] do
    Y_Opi := [];
    Xi := [];
    if elleqoo then
      Append(~Xi, 1);
      Append(~Y_Opi, Op!(lambda*lifts[i]*lambda^(-1)));
    else
      if inNormSupport then
        // Work hard
        for j in [1..#lambdas] do
          for c in [1..#lambdas] do
            if lambdas[j]*lifts[i]*(lambdas[c])^(-1) in Op then
              Append(~Xi, c);
              Append(~Y_Opi, Op!(lambdas[j]*lifts[i]*lambdas[c]^(-1)));
              break c;
            end if;
          end for;
        end for;
      else
        for j in [1..numP1] do
          iotadelta := iotaell(ellcosets[j]*lifts[i]);
          _, v := P1ellrep(iotadelta[2], false, false);
          c := Index(P1ell, v);
          Append(~Xi, c);
          Append(~Y_Opi, Op!(levelmults[j]*lambda*ellcosets[j]*lifts[i]*
                             (levelmults[c]*lambda*ellcosets[c])^(-1)));
        end for;
      end if;
    end if;
    Append(~X, Xi);
    Append(~Y_Op, Y_Opi);
  end for;

// GetMemoryUsage(); MemProfile();

  Y_U := [];
  vprintf ModFrmHil: "Reducing %4o units of Gamma ......................... ", n*numP1;
  vtime ModFrmHil:
  for i in [1..n] do
    G := [];

    for j in [1..numP1] do
      y := CompleteRelationFromUnit(Gammap, Y_Op[i][j], RPAsp, RPAspinv : IsTrivialCoefficientModule := IsLevelOne);
      if not IsLevelOne then
        y := y*Zp[X[i][j]];
      end if;
      Append(~G, y);
    end for;
    Append(~Y_U, G);
  end for;

// GetMemoryUsage(); MemProfile();

  vprintf ModFrmHil: "Computing H1 (coinduced) ............................. ";
  vtime ModFrmHil:
  Htilde, mH := InducedH1(Gamma, Gammap, N, cosets, cosetsp, RPAs, RPAsinv, RPAsp, RPAspinv);
  
  if #Htilde eq 0 then
    return [];
  else
    M := HorizontalJoin([ &+[ Y_U[i][j] : j in [1..numP1]] : i in [1..n] ]);
    MH := Matrix(Htilde)*M;
    MM := Matrix([mH(MH[i]) : i in [1..#Htilde]]);
    return MM;
  end if;
end function;
