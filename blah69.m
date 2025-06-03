load "config.m";

SetVerbose("ModFrmHil", 3);

F := QuadraticField(5);
ZF := Integers(F);

epses := UnitsGenerators(F : exclude_torsion:=false);
vs := InfinitePlaces(F)[1];
k := [1,2];

for N in [N : N in IdealsUpTo(1000, F) | not IsZero(N)] do
  H := HeckeCharacterGroup(N, [1,2]);
  for chi in Elements(H) do
    if IsCompatibleWeight(chi, k) then
      chi_ct +:= 1;
      print "---------------------------------";
      print "N", IdealOneLine(N);
      print "chi", chi, "order", Order(chi);
      print HeckeCharLabel(chi);
    end if;
  end for;
end for;
