// this test takes too long, so I've disabled it, but
// it computes a weight [2, 3] space with nebentypus
// and checks that the square of the space lies within
// the weight [4, 6] forms.

F := QuadraticField(2);
ZF := Integers(F);

M := GradedRingOfHMFs(F, 500);
N := ideal<ZF | 119, -ZF.2 - 45>;
k := [1, 2];

H := HeckeCharacterGroup(N, [1, 2]);
chi := H.2;
assert HeckeCharLabel(chi_12) eq "-2.0.1_119.2_2u1.0u1.2u";

Mk := HMFSpace(M, N, [2, 3], chi);
Sk := CuspFormBasis(Mk);

M46 := HMFSpace(M, N, [4, 6]);
S46 := CuspFormBasis(M46);

#LinearDependence(Sk) eq 0;
#LinearDependence(S46) eq 0;
f := Sk[1];
#Intersection([f^2], S46) eq 1;
