/******************************************************************************
* This file has tests for ModFrmHil/level.m::HeckeMatrix2, a modification of 
* !Geometry/ModFrmHil/level.m::HeckeMatrix.
*
* The tests ensure that the modifications do not break any existing behavior.
*******************************************************************************/

MAX_PP_NORM := 50;
NUM_TRIALS := 2;
OPTIONAL_TESTS := false;

// set to true if you want to print times
PRINT_TIMES := true;

import "./ModFrmHil/ideal_datum.m" : induced_module_mtrxs_of_gens;

forward RightPermutationActions;


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
  return RPAs;
end function;

// compares HeckeMatrix(Gamma, N, pp) to HeckeMatrix2(Gamma, N, pp)
// and optionally returns their respective durations

// Note that the first iteration of HeckeMatrix in each test() may be
// very slow if a presentation for Gamma needs to be computed. Once this
// is computed it is cached and subsequent runs can be meaningfully
// compared.
procedure compare_and_time(Gamma, N, pp, k : UseAtkinLehner:=false)
  F := BaseField(QuaternionAlgebra(Gamma));

  t := Cputime();
  M2 := HeckeMatrix2(Gamma, N, pp, k, 1 : UseAtkinLehner:=UseAtkinLehner);
  hm2_time := Cputime(t);

  t := Cputime();
  M := HeckeMatrix(Gamma, N, pp : UseAtkinLehner:=UseAtkinLehner);
  hm_time := Cputime(t);

  assert M eq M2;

  if PRINT_TIMES then
    print "HeckeMatrix time: ", hm_time;
    print "HeckeMatrix2 time: ", hm2_time;
    print "\n";
  end if;
end procedure;


// the default quaternion algebra is the one ramified at all but one infinite place
// this will throw an error if B over a field of even degree.
procedure test(F, N : B:=0)
  // F::FldNum
  // N::RngOrdIdl

  if PRINT_TIMES then
    print "-----------------------------------------------\n";
  end if;

  ZF := Integers(F);

  if B cmpeq 0 then
    B := QuaternionAlgebra(1*ZF, [InfinitePlaces(F)[i] : i in [2 .. Degree(F)]]);
  end if;

  assert BaseField(B) eq F;
  O := MaximalOrder(B);
  D := Discriminant(B);
  // the level and discriminant are assumed to be coprime
  assert IsCoprime(N, D);
  Gamma := FuchsianGroup(O);

  X := cIdealDatum(Gamma, N);
  // because HeckeMatrix only works in parallel weight, this is the setting
  // in which we test
  k := [2 : _ in [1 .. Degree(F)]];

  assert [<a, b> : a->b in induced_module_mtrxs_of_gens(X, k)] eq 
              [<a, b> : a->b in RightPermutationActions(X)];

  // Hecke matrix at infinity
  compare_and_time(Gamma, N, "Infinity", k);

  coprime_pps := [pp : pp in PrimesUpTo(MAX_PP_NORM, F) | GCD(pp, N * D) eq 1*ZF];
  level_pps := [fac[1] : fac in Factorization(N)];
  disc_pps := [fac[1] : fac in Factorization(Discriminant(B))];

  for _ in [1 .. NUM_TRIALS] do
    pp := Random(coprime_pps);
    compare_and_time(Gamma, N, pp, k);

    if #level_pps ne 0 then
      pp := Random(level_pps);
      compare_and_time(Gamma, N, pp, k);
      compare_and_time(Gamma, N, pp, k : UseAtkinLehner:=true);
      Exclude(~level_pps, pp);
    end if;

    if #disc_pps ne 0 then
      pp := Random(disc_pps);
      compare_and_time(Gamma, N, pp, k : UseAtkinLehner:=true);
      Exclude(~disc_pps, pp);
    end if;
  end for;
end procedure;

R<x> := PolynomialRing(Rationals());

// h+(F) = 1, D = 1, N = 4
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

test(F, 4*ZF);

// h+(F) = 1, D = 1, N = (7, a)^2 * (2)
pps := PrimesUpTo(15, F);
test(F, pps[1]^2 * pps[2]);

// h+(F) = 1, D = (2) * (13, alpha), N = (7, a)
B := QuaternionAlgebra(2 * Factorization(13*ZF)[1][1], [InfinitePlaces(F)[i] : i in [2 .. Degree(F)]]);
test(F, Factorization(7*ZF)[1][1] : B:=B);

// h+(F) = 3, D = 1
F := NumberField(x^3-x^2-9*x+8);
ZF := Integers(F);
pps := PrimesUpTo(15, F);
test(F, pps[1] * pps[2]);

if OPTIONAL_TESTS then
  // h+(F) = 1, D = 1, N = 1
  F := NumberField(x^3-x^2-10*x+8);
  ZF := Integers(F);
  // the dimension of H1 is 4, see Greenberg-Voight
  test(F, 1*ZF);

  // h+(F) = 1, D = 1, N = 3
  test(F, 3*ZF);

  // h+(F) = 2, D = 1
  F := NumberField(x^3-x^2-9*x+10);
  ZF := Integers(F);
  test(F, 2*ZF);
end if;

