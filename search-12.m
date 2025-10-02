load "config.m";

F := QuadraticField(2);
ZF := Integers(F);
M := GradedRingOfHMFs(F, 777);
k := [1, 2];

ideals := [N : N in IdealsUpTo(600, F) | Norm(N) eq StringToInteger(N_norm)];

for N in ideals do
  H := HeckeCharacterGroup(N, [1,2]);
  for chi in Elements(H) do
    if IsCompatibleWeight(chi, k) and Order(chi) eq 2 then
      chi_label := HeckeCharLabel(chi);
      print "--------------------- Trying:", IdealOneLine(N), chi, chi_label;
      cond_fin, cond_inf := Conductor(chi);
      Mk := HMFSpace(M, N, k, chi);
      try
        Sk := CuspFormBasis(Mk : StableOnly:=true);
        if #Sk gt 0 then
          print "******************** Found a [1,2] form:", IdealOneLine(N), chi, chi_label;
        end if;
      catch e;
        print "!!!!!!!!!!!!!!!!!!!!!! error:", IdealOneLine(N), chi, chi_label;
        print e;
      end try;
    end if;
  end for;
end for;
