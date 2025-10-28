load "config.m";

F := QuadraticField(2);
ZF := Integers(F);
M := GradedRingOfHMFs(F, 500);
k := [1, 2];

ideals := [N : N in IdealsUpTo(1000, F) | Norm(N) ge 195];
  print Norm(N);
  H := HeckeCharacterGroup(N, [1,2]);
  chis := [chi : chi in Elements(H) | Order(chi) eq 2 and IsCompatibleWeight(chi, k)];
  for chi in chis do
    Mk := HMFSpace(M, N, k, chi);
    num_gcs := #PossibleGrossenchars(Mk);
    if num_gcs gt 0 then
      print "Found:", LMFDBLabel(N), HeckeCharLabel(chi), chi, num_gcs;
      // print Factorization(N);
      // Sk := CuspFormBasis(Mk : StableOnly:=true);
      // print "#Sk", #Sk;
    end if;
  end for;
end for;
