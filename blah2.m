// IF YOU SEE THIS AGAIN, THIS IS A BUG IN THE DIHEDRAL CODE AND NEEDS
// TO BE FIXED AT SOME POINT!!!!!!!!

load "config.m";
F := QuadraticField(5);
ZF := Integers(F);
N := 3^2 * Factorization(5*ZF)[1][1]^2;
print N, Norm(N);
k := [2,2];
GRng := GradedRingOfHMFs(F, 150);
H := HeckeCharacterGroup(N, [1,2]);
chi := H.0;
Mk := HMFSpace(GRng, N, k, chi);
// psis := PossibleGrossenchars(Mk);

Dk := DihedralBasis(Mk);


