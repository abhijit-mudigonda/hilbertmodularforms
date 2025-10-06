// find (N, chi) tuples supporting dihedral [1,2] forms

load "config.m";

// SetVerbose("ModFrmHil", 3);
// SetVerbose("HilbertModularForms", 3);

F := QuadraticField(2);
ZF := Integers(F);
M := GradedRingOfHMFs(F, 300);
k := [3, 3];

/*
N := 5*ZF;
H := HeckeCharacterGroup(N, [1,2]);
chi := H.1;
assert HeckeCharLabel(chi) eq "-2.0.1_25.1_4u1u1.2u";

Mk := HMFSpace(M, N, k, chi);
Sk := CuspFormBasis(Mk);
Sk_hs := HeckeStabilityCuspBasis(Mk : prove:=false, stable_only:=true);
*/

ideals := [N : N in IdealsUpTo(500, F) | N eq 1*ZF or Max([tup[2] : tup in Factorization(N)]) eq 1];
for N in ideals do
  print Norm(N), IdealOneLine(N);
  H := HeckeCharacterGroup(N, [1,2]);
  chis := [chi : chi in Elements(H) | IsCompatibleWeight(chi, k) and Order(chi) gt 2];
  for chi in chis do
    Mk := HMFSpace(M, N, k, chi);
    if #HeckeStabilityCuspBasis(Mk : prove:=false, stable_only:=true) gt 0 then
      print "hi!!!", chi, Order(chi), HeckeCharLabel(chi);
    end if;
  end for;
end for;
