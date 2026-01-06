load "config.m";

MAX_PRIME_NORM := 500;

F := QuadraticField(2);
ZF := Integers(F);
M := GradedRingOfHMFs(F, 800);
k := [1, 2];
F_aut := Automorphisms(F)[2];

full_labels := [
  "-2.0.1_119.2_2u1.0u1.2u",
  "-2.0.1_196.1_2u0.1u1.2u",
  "-2.0.1_224.1_2u1.0.0.1u1.2u",
  "-2.0.1_252.1_2u1.1.0u1.2u",
  "-2.0.1_272.2_2u1.1.0u1.2u",
  "-2.0.1_272.1_2u0.0.1u1.2u"
  ];

/*
for full_label in full_labels do
  print "*************************************************", full_label;
  process_full_label(full_label);
end for;
*/

// full_label := "-2.0.1_512.1_2u0.1.1u1.2u";
full_label := "-2.0.1_119.2_2u1.0u1.2u";

chi := FullChiLabelToHeckeChar(full_label);
N := Modulus(chi);
Mk := HMFSpace(M, N, k, chi);

K_rels := [K_rel : K_rel in QuadraticExtensionsWithConductor(N, [1,2]) | IsTotallyComplex(K_rel)];
psis := [* *];

for K_rel in K_rels do
  psis cat:= [* psi : psi in PossibleGrossencharsOfRelQuadExt(K_rel, N, k, chi) *];
end for;

psis := PossibleGrossenchars(Mk);
assert #psis eq 2;
psi := psis[1];
X := Parent(psi);

K := X`BaseField;
ZK := Integers(K);

custom_wt := HeckeCharWeightFromWeight(K, F, [0, 2]);

for p in PrimesUpTo(MAX_PRIME_NORM) do
  assert Type(chi(p)) eq RngIntElt;
end for;

twisted_apps := AssociativeArray();

for pp in PrimesUpTo(MAX_PRIME_NORM, F) do
  fact := Factorization(ZK !! pp);
  g := #fact;
  d := InertiaDegree(pp);
  if g eq 2 then
    twisted_apps[pp] := StrongAdd(
      [* 
        Evaluate(psi, fact[1][1] : custom_weight:=custom_wt),
        Evaluate(psi, fact[2][1] : custom_weight:=custom_wt) 
      *]);
  elif fact[1][2] ne 1 then
    twisted_apps[pp] := Evaluate(psi, fact[1][1] : custom_weight:=custom_wt);
  else
    twisted_apps[pp] := 0;
  end if;
end for;

split_primes := [
  <Factorization(ZF!!p)[1][1], Factorization(ZF!!p)[2][1]>
  : p in PrimesUpTo(MAX_PRIME_NORM) | #Factorization(ZF!!p) eq 2
  ];

aps := AssociativeArray();

for pp_tup in split_primes do
  print "---", Norm(pp_tup[1]);
  prod_apps := StrongMultiply([* twisted_apps[pp_tup[1]], twisted_apps[pp_tup[2]] *]);
  b, x := IsCoercible(Rationals(), prod_apps);
  if b then
    print x;
  else
    print b;
    b1, app1_squared := IsCoercible(F, twisted_apps[pp_tup[1]]^2);
    b2, app2_squared := IsCoercible(F, twisted_apps[pp_tup[2]]^2);

    print b1 and b2;
    if b1 and b2 then
      print F_aut(app1_squared) eq app2_squared;
    end if;
  end if;
end for;

psi := psis[1];
f := ThetaSeries(Mk, psi);

M24 := HMFSpace(M, N, [2,4]);
S24 := CuspFormBasis(M24);
