load "config.m";

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

M := GradedRingOfHMFs(F, 150);

k := [1,1,2];
ideals := [N : N in IdealsUpTo(3000, F) | not IsZero(N)];

for N in ideals do
  H := HeckeCharacterGroup(N, [1,2,3]);
  for chi in Elements(H) do
    if IsCompatibleWeight(chi, k) and Order(chi) eq 2 then
      print "--------------------- Trying:", IdealOneLine(N), chi, HeckeCharLabel(chi);
      Mk := HMFSpace(M, N, k, chi);
      try
        Sk := HeckeStabilityCuspBasis(Mk : prove:=false, stable_only:=true);
        if #Sk gt 0 then
          print "******************** Found a [1,1,2] form:", IdealOneLine(N), chi, HeckeCharLabel(chi);
        end if;
      catch e;
        print e;
      end try;
    end if;
  end for;
end for;
