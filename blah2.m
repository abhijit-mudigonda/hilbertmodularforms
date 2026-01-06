load "config.m";

SetVerbose("ModFrmHil", 3);
SetVerbose("HilbertModularForms", 3);

F := QuadraticField(2);
ZF := Integers(F);
k := [2,3];

chi_label := "-2.0.1_391.1_2u1.1u1.2u";
qq := ideal<ZF | 17, ZF.2 + 6>;

chi := FullChiLabelToHeckeChar(chi_label);
assert Order(chi) eq 2;
assert IsCompatibleWeight(chi, k);
N := Modulus(chi);

M := GradedRingOfHMFs(F, 200);
M23 := HMFSpace(M, N, k, chi);
S23 := CuspFormBasis(M23);

print "done S23, now S46!";
M46 := HMFSpace(M, N, [4,6]);
S46 := CuspFormBasis(M46);

S23_squared := Basis(&cat[[f * g : g in S23] : f in S23]);

print "dim S23", #S23;
print "dim S23_squared", #S23_squared;
print "dim S46", #S46;
print "intersection", #Intersection(S23_squared, S46);

