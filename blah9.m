load "config.m";

F := QuadraticField(5);
ZF := Integers(F);
M := GradedRingOfHMFs(F, 444);
k := [2, 2];
N := ideal<ZF | 38, 4*ZF.2 + 16>;
H := HeckeCharacterGroup(N, [1,2]);
chi := H.1^2;
assert Order(chi) eq 3;
Mk := HMFSpace(M, N, k, chi);
Sk := CuspFormBasis(Mk);


/*
for N in IdealsUpTo(100, F) do
  print Norm(N), IdealOneLine(N);
  H := HeckeCharacterGroup(N, [1,2]);
  chis := [chi : chi in Elements(H) | Order(chi) gt 2 and IsCompatibleWeight(chi, k)];
  for chi in chis do
    Mk := HMFSpace(M, N, k, chi);
    print "chi", chi;
    Sk := CuspFormBasis(Mk);
    if #Sk gt 0 then
      print "found one!", chi, HeckeCharLabel(chi), #Sk;
    end if;
  end for;
end for;
*/
