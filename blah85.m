load "config.m";

F := QuadraticField(5);
k := [1,9];
M := GradedRingOfHMFs(F, 100);
for N in PrimesUpTo(500, F) do
  print IdealOneLine(N);
  H := HeckeCharacterGroup(N, [1,2]);
  for chi in Elements(H) do
    if IsCompatibleWeight(chi, k) then
      Mk := HMFSpace(M, N, k, chi);
      Dk := DihedralBasis(Mk);
      if #Dk ne 0 then
        print "Found one!";
        print N, chi;
      end if;
    end if;
  end for;
end for;
