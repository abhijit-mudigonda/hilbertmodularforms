/*************************************************************
 * Computing a 2-dimensional space of nonparitious forms
 * of level norm 7 and weight (4, 3). 
 * This reproduces an example first computed by
 * Dembélé-Loeffler-Pacetti (https://arxiv.org/pdf/1612.06625)
**************************************************************/

t := Cputime();

load "config.m";

BOUND := 200;

// specify the field, level, weight, and nebentypus
F := QuadraticField(2);
ZF := Integers(F);
k := [4, 3];
N := Factorization(7*ZF)[1][1];
H := HeckeCharacterGroup(N, [1,2]);
chi := H.1;

// compute an basis of q-expansions for the corresponding space
// the basis will be over F'(chi) = F
M := GradedRingOfHMFs(F, BOUND);
Mk := HMFSpace(M, N, k, chi);
Sk := CuspFormBasis(Mk);

// the dimension of the cusp space is 2
assert #Sk eq 2;

// diagonalizing the action of Hecke operators on the space
// to produce an eigenbasis
eigs := Eigenbasis(Mk, Sk : P:=10);

print Cputime(t);

f := eigs[1];
bb := 1*ZF;

coeffs := AssociativeArray();
dinv := M`IdealToRep[bb][1*ZF];
f := f / Coefficient(f, bb, dinv);

for pp in PrimesUpTo(BOUND, F) do
  pidinv := M`IdealToRep[bb][pp];
  coeffs[pp] := Coefficient(f, bb, pidinv);
  print Norm(pp), "|", pidinv / dinv, "|", Coefficient(f, bb, pidinv);
end for;
