// check for the CM (1,2) form at level norm 391
// computes the space, checks that we see the CM forms, 
// checks that the squares are in (2,4), and also does some 
// oldform newform checks since that's where the issues were.

load "config.m";

SetVerbose("HilbertModularForms", 3);
// SetVerbose("ModFrmHil", 3);

F := QuadraticField(2);
M := GradedRingOfHMFs(F, 200);
k := [1, 2];

chi_label := "-2.0.1_391.1_2u1.0u1.2u";
chi := FullChiLabelToHeckeChar(chi_label);
N := Modulus(chi);

Mk := HMFSpace(M, N, k, chi);
// Sk := CuspFormBasis(Mk : StableOnly:=true);

Dk := DihedralBasis(Mk);

M24 := HMFSpace(M, N, [2,4]);
// S24 := CuspFormBasis(M24);
Dk_squared := [Dk[1]^2, Dk[2]^2, Dk[1] * Dk[2]];

H := Parent(chi);
chi_eis := H.1 * H.2;
M11 := HMFSpace(M, N, [1,1], chi_eis);
h := EisensteinBasis(M11)[1];
chi_23 := H.2;

M23 := HMFSpace(M, N, [2,3], chi_23);
S23 := CuspFormBasis(M23);
assert #LinearDependence(S23) eq 0;
#LinearDependence(S23 cat [Dk[1] * h]);

/*
new_basis := NewCuspFormBasis(M23);
old_basis := OldCuspFormBasis(M23);

assert #LinearDependence(new_basis cat old_basis) eq 0;

[#Intersection(old_basis, [HeckeOperator(f, pp) : f in old_basis]) : pp in PrimesUpTo(10, F)];
*/
