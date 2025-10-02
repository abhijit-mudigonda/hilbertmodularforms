load "config.m";

F := QuadraticField(2);
ZF := Integers(F);
M := GradedRingOfHMFs(F, 10);
k := [1, 2];

for N in IdealsUpTo(1000, F) do
  print Norm(N);
  H := HeckeCharacterGroup(N, [1,2]);
  chis := [chi : chi in Elements(H) | Order(chi) eq 2 and IsCompatibleWeight(chi, k)];
  for chi in chis do
    Mk := HMFSpace(M, N, k, chi);
    if #PossibleGrossenchars(Mk) gt 0 then
      print "Found:", LMFDBLabel(N), HeckeCharLabel(chi), chi;
    end if;
  end for;
end for;
