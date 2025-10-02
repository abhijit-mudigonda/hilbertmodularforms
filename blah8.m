// find (N, chi) tuples supporting dihedral [1,2] forms

load "config.m";

F := QuadraticField(2);
ZF := Integers(F);
M := GradedRingOfHMFs(F, 300);
k := [1, 2];

for N in IdealsUpTo(500, F) do
  print Norm(N), IdealOneLine(N);
  H := HeckeCharacterGroup(N, [1,2]);
  chis := [chi : chi in Elements(H) | IsCompatibleWeight(chi, k)];
  for chi in chis do
    Mk := HMFSpace(M, N, k, chi);
    if #PossibleGrossenchars(Mk) gt 0 then
      print "hi!!!", chi, Order(chi), HeckeCharLabel(chi);
    end if;
  end for;
end for;

