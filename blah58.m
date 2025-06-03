load "config.m";

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

M := GradedRingOfHMFs(F, 150);

k := [2,2,3];
ideals := [I : I in IdealsUpTo(600, F)];

for N in ideals do
  H := HeckeCharacterGroup(N, [1,2,3]);
  chis := [chi : chi in Elements(H) | Order(chi) eq 2];
  if &or[IsCompatibleWeight(chi, [1,1,1]) : chi in chis] and &or[IsCompatibleWeight(chi, [1,1,2]) : chi in chis] then
    print IdealOneLine(N);
  end if;
  /*
  for chi in Elements(H) do
    if IsCompatibleWeight(chi, k) and Order(chi) eq 2 then
      print "--------------------- Trying:", IdealOneLine(N), chi, HeckeCharLabel(chi);
      Mk := HMFSpace(M, N, k, chi);
      Sk := HeckeStabilityCuspBasis(Mk : prove:=false, stable_only:=true);
      if #Sk gt 0 then
        print "******************** Found a [1,1,2] form:", IdealOneLine(N), chi, HeckeCharLabel(chi);
      end if;
    end if;
  end for;
  */
end for;
