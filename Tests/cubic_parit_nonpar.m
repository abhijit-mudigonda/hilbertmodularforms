R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

prec := 300;
M := GradedRingOfHMFs(F, prec);
N := ideal<ZF | 6*ZF.1 - 3*ZF.2 - 2*ZF.3>;
k := [1,1,3];
chi := HeckeCharacterGroup(N, [1,2,3]).1;

Mk := HMFSpace(M, N, k, chi);
Sk := HeckeStabilityCuspBasis(Mk : prove:=false);
Dk := DihedralBasis(Mk);

assert #Sk eq 1;
assert #Dk eq 1;

// the unique cusp form in this space is dihedral
assert #LinearDependence([Sk[1], Dk[1]]) eq 1;

