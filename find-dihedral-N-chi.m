load "config.m";

F := QuadraticField(2);
ZF := Integers(F);
M := GradedRingOfHMFs(F, 10);
k := [1, 2];

for N in IdealsUpTo(1000, F) do
  H := HeckeCharacterGroup(N, [1,2]);
  chis := [chi : chi in Elements(H) | Order(chi) eq 2 and IsCompatibleWeight(chi, k)];
  for chi in chis do
    Mk := HMFSpace(M, N, k, chi);
    if #PossibleGrossenchars(Mk) gt 0 then
      flag := false;
      flag_quad := false;
      for psi in Elements(H) do
        if IsCompatibleWeight(psi, [1, 1]) and IsGamma1EisensteinWeight(psi, 1) then
          chi_eis := psi;
          flag := true;
          if Order(psi) eq 2 then
            flag_quad := true;
          end if;
        end if;
      end for;
      print LMFDBLabel(N), HeckeCharLabel(chi), chi, flag, flag_quad;
    end if;
  end for;
end for;
