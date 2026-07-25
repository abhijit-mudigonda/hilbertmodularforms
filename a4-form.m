// This is code for producing a parallel weight 1 A4-form over F = Q(\sqrt{5})
// and also for producing a degree 8 field field K/Q such that the splitting field 
// of K/F is the A4 extension through which the projective representation factors.
// It also produces some data that could be useful for computing the actual group
// through which the extension factors. 

load "config.m";
// SetVerbose("HilbertModularForms", 3);

// MAX_NORM := 150;
F := QuadraticField(5);
ZF := Integers(F);

function orderer(x)
  if x eq 4 then
    return 1;
  elif x eq 0 then
    return 2;
  elif x eq 1 then
    return 3;
  elif x eq 2 then
    return 4;
  elif x eq (3 + F.1)/2 then
    return 5;
  elif x eq 3 then
    return 6;
  else
    print "busted!";
    assert 0 eq 1;
  end if;
end function;

PREC := 1000;
M := GradedRingOfHMFs(F, PREC);
N := 18*Factorization(5*ZF)[1][1];
k := [1,1];


H := HeckeCharacterGroup(N, [1,2]);
chis := [chi : chi in Elements(H) | IsCompatibleWeight(chi, k) and (Order(chi) eq 6)];

primes := PrimesUpTo(PREC, F : coprime_to:=N);

M22 := HMFSpace(M, N, [2,2]);
S22 := CuspFormBasis(M22);

print "Computed a basis for the cusp forms of weight 2!";


chi := H.2 * H.4;
Mk := HMFSpace(M, N, k, chi);
Sk := CuspFormBasis(Mk : stable_only:=true);
Ek := EisensteinBasis(Mk);
eigs := Eigenbasis(Mk, Sk : P:=20);
eigs := [DivideByFirstNonzeroIdlCoeff(f) : f in eigs];

f := eigs[1];
assert #Intersection([f], Ek) eq 0;
assert not ProbabilisticDihedralTest(f);

f := ChangeCoefficientRing(f, CyclotomicField(60));
E := CyclotomicField(12);
f := ChangeCoefficientRing(f, E);


// the a_pp as well as the order of Frob_pp in the projective rep
A := [<Norm(pp), IdealOneLine(pp), Coefficient(f, pp), orderer(Integer\
s()!(Coefficient(f, pp)^2 / chi(pp)))> : pp in primes];
// the conjugacy classes in the (non-projective) Galois image
B := Multiset([<Coefficient(f, pp), E!chi(pp)> : pp in primes]);

root_dict := AssociativeArray();
for a in [0 .. 11] do
  root_dict[E.1^a] := a;
end for;

S<y> := PolynomialRing(E);
// stores the eigenvalues of the various conjugacy classes
C := Multiset([[root[1] : root in Roots(y^2 - Coefficient(f, pp)*y + E!chi(pp))] : pp in primes]);
// shows the orders of the elements - should show elements of order 2, 3, 4, 6, and 12 at precision 800
D := Multiset([LCM([Integers()!(12 / GCD(12, root_dict[root[1]])) : root in Roots(y^2 - Coefficient(f, pp)*y + E!chi(pp))]) : pp in primes]);


// The number field corresponding to a degree 8 subfield
// of the A4 extension through which the projective rep factors
//
// This is a degree 8 subfield fixed by a 3-cycle
R<x> := PolynomialRing(Rationals());
K<a> := NumberField(x^8 - 3*x^7 + 4*x^6 - 7*x^5 + 6*x^4 + 5*x^3 + 5*x^2 - 30*x + 20);
F := QuadraticField(5);

_ := IsSubfield(F, K);
K_rel := RelativeField(F, K);
assert Degree(K_rel) eq 4;

ZK_rel := Integers(K_rel);

// primes of F splitting into two primes with inertia
// degrees (3, 1), two primes with inertia degrees (2, 2)
// and four primes, respectively
three_ones := [];
two_twos := [];
fours := [];

for pp in primes do
  facts := Factorization(ZK_rel!!pp);
  if #facts eq 4 then
    Append(~fours, pp);
  elif #facts eq 2 then
    if Norm(facts[1][1]) eq Norm(facts[2][1]) then
      Append(~two_twos, pp);
    else
      Append(~three_ones, pp);
    end if;
  end if;
end for;

// this is checking that the projective Galois rep
// and the splitting field of the degree 8 extension
// actually agree
assert [Norm(pp) : pp in three_ones] eq [a[1] : a in A | a[4] eq 3];
assert [Norm(pp) : pp in two_twos] eq [a[1] : a in A | a[4] eq 2];
assert [Norm(pp) : pp in fours] eq [a[1] : a in A | a[4] eq 1];

/*

// this is code to search for the A4 form, not relevant now but might help later
primes := PrimesUpTo(MAX_NORM, F : coprime_to:=N);

Mks := AssociativeArray();
Sks := AssociativeArray();

print "number of chis is", #chis;
for chi in chis_2 do
  print "--------------------------------------";
  print "processing chi", chi;
  Mk := HMFSpace(M, N, k, chi);
  Sk := CuspFormBasis(Mk : stable_only:=true);
  Ek := EisensteinBasis(Mk);
  Sks[chi] := Sk;
  Mks[chi] := Mk;

  if #Sk - #Intersection(Sk, Ek) gt 0 then
    if #Sk gt 1 then
      B := Eigenbasis(Mk, Sk : P:=17);
    else
      B := Sk;
    end if;

    B := [DivideByFirstNonzeroIdlCoeff(f) : f in B];
    for f in B do
      if (not ProbabilisticDihedralTest(f)) and (#Intersection([f], Ek) eq 0) then
        print "Found a form that's neither dihedral nor Eisenstein";
        orders := [orderer(Integers()!(Coefficient(f, pp)^2 / chi(pp))) : pp in primes];
        print orders;
        if not (6 in orders) then
          print "This one might be non-dihedral!", chi;
        end if;
        print [<Norm(primes[i]), IdealOneLine(primes[i]), orders[i]> : i in [1..#primes]];
      end if;
    end for;
  end if;
end for;
*/
