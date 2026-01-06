load "config.m";

F := QuadraticField(5);
M := GradedRingOfHMFs(F,
k := [3,5];
for N in [I : I in IdealsUpTo(100, F) | not IsZero(N)] do
  for chi in HeckeCharacterGroup(N, [1,2]) do
    if IsCompatibleWeight(chi, k) then
      Mk := HMFSpace(M, N, k, chi);
      Sk := CuspFormBasis(Mk);
      if #Sk gt 0 then
        print "found one!", IdealOneLine(N), chi, HeckeCharLabel(chi);
      end if;
    end if;
  end for;
end for;
