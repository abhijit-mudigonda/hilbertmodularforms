load "config.m";

/* 
// https://www.lmfdb.org/ModularForm/GL2/Q/holomorphic/256/2/a/e/
K := QuadraticField(-2);
ZK := Integers(K);
N := Factorization(2*ZK)[1][1]^5;
*/

/* idr what this is
K := QuadraticField(-2);
ZK := Integers(K);
pp := Factorization(3*ZK)[1][1];
assert Norm(pp) eq 3;
N := pp;
*/

// https://www.lmfdb.org/ModularForm/GL2/Q/holomorphic/484/2/a/d/
K := QuadraticField(-11);
ZK := Integers(K);
N := 2*ZK * Factorization(11*ZK)[1][1];

D := DirichletGroup(N);
ps := [p : p in PrimesUpTo(70) | #Factorization(p*ZK) gt 1];
aut := Automorphisms(K)[2];
psis := Elements(D);

assert NarrowClassNumber(K) eq 1;

pis := [];
for p in ps do
  _, pi := IsNarrowlyPrincipal(Factorization(p*ZK)[1][1]);
  Append(~pis, pi);
end for;


function aps_from_psi(psi)
  L := Compositum(CyclotomicField(Order(psi)), K);
  a_ps := [];
  for pi in pis do
    a_p := L!pi * L!psi(pi) + L!aut(pi) * L!psi(aut(pi));
    Append(~a_ps, <Norm(pi), a_p>);
  end for;
  return a_ps;
end function;

function eval_hc(psi, pp)
  L := Compositum(CyclotomicField(Order(psi)), K);
  b, pi := IsNarrowlyPrincipal(pp);
  assert b;
  return pi * psi(pi);
end function;

psis_aps := [* *];
for psi in psis do
  Append(~psis_aps, <psi, aps_from_psi(psi)>);
end for;

/* 
// https://www.lmfdb.org/ModularForm/GL2/Q/holomorphic/256/2/a/e/

psi := psis[6];
aps := aps_from_psi(psi);
aps;

ap_ress := [];
for ap_tup in aps do
  p, ap := Explode(ap_tup);
  if p mod 4 eq 1 then
    Append(~ap_ress, <p, ap>);
  elif p mod 4 eq 3 then
    Append(~ap_ress, <p, ap^2>);
  end if;
end for;

ap_ress;

*/
    

