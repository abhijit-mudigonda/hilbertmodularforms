load "config.m";

F := QuadraticField(5);
M := GradedRingOfHMFs(F, 50);

k := [4,4];

for N in IdealsUpTo(20, F) do
  Mk := HMFSpace(M, N, k);
  Sk := #CuspFormBasis(Mk);
  print IdealOneLine(N), "#Sk", Sk;
end for;

