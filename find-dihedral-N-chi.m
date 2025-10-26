load "config.m";

F := QuadraticField(2);
ZF := Integers(F);
M := GradedRingOfHMFs(F, 500);
k := [1, 2];

ideals := [N : N in IdealsUpTo(1000, F) | Norm(N) ge 195];
for N in ideals do
  H := HeckeCharacterGroup(N, [1,2]);
  chis := [chi : chi in Elements(H) | Order(chi) eq 2 and IsCompatibleWeight(chi, k)];
  for chi in chis do
    Mk := HMFSpace(M, N, k, chi);
    if #PossibleGrossenchars(Mk) gt 0 then
      print IdealOneLine(N), HeckeCharLabel(chi), chi;
      Sk := CuspFormBasis(Mk : StableOnly:=true);
      print "possible grossenchars", #PossibleGrossenchars(Mk);
      print "#Sk", #Sk;
    end if;
  end for;
end for;
