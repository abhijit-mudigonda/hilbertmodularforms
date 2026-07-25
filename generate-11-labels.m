load "config.m";

LOW_N_NORM := 1620;
HIGH_N_NORM := 4000;

F := QuadraticField(5);
ZF := Integers(F);
k := [1, 1];

out := Open("11-labels-1620-4000.txt", "w");
for N in IdealsUpTo(HIGH_N_NORM, F) do
  if Norm(N) ge LOW_N_NORM then
    H := HeckeCharacterGroup(N, [1,2]);
    for chi in Elements(H) do
      if IsCompatibleWeight(chi, k) then
        fprintf out, "%o\n", HeckeCharLabel(chi);
      end if;
    end for;
  end if;
end for;
delete out;
