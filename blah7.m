// check for the CM (1,2) form at level norm 441
// computes the space, checks that we see the CM forms, 
// checks that the squares are in (2,4), and also does some 
// oldform newform checks since that's where the issues were.

load "config.m";

SetVerbose("HilbertModularForms", 3);
SetVerbose("ModFrmHil", 3);

F := QuadraticField(5);
ZF := Integers(F);
M := GradedRingOfHMFs(F, 888);
k := [1, 2];

chi_label := "-5.0.1_576.1_2u0.1.0.1u1.2u"
chi := FullChiLabelToHeckeChar(chi_label);
N := Modulus(chi);

Mk := HMFSpace(M, N, k, chi);
Sk := CuspFormBasis(Mk : StableOnly:=true);

print "dim of Sk is", #Sk;

pp := PrimesUpTo(20, F)[2];
V := HeckeStableSubspace(Sk, pp : use_fourier:=true);

M24 := HMFSpace(M, N, [2,4]);
S24 := CuspFormBasis(M24);
assert #LinearDependence(S24) eq 0;

Sk_squared := Basis(&cat[[f * g : g in Sk] : f in Sk]);
print "dim of Sk_squared is", #Sk_squared;
print "dim of Sk_squared in S24 is", #Intersection(Sk_squared, S24);

/*
new_basis := NewCuspFormBasis(M23);
old_basis := OldCuspFormBasis(M23);

assert #LinearDependence(new_basis cat old_basis) eq 0;

[#Intersection(old_basis, [HeckeOperator(f, pp) : f in old_basis]) : pp in PrimesUpTo(10, F)];
*/

// hm := HeckeMatrix(S23, Factorization(7*ZF)[1][1] : use_coeff:=true);
