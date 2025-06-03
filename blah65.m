load "config.m";

SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

N := ideal<ZF | [1169, 0, 0], [291, 0, -1]>;

epses := UnitsGenerators(F : exclude_torsion:=false);
vs := InfinitePlaces(F)[1 .. 2];
k := [1,1,2];

f := Open("N_chis_all.txt", 'w');

for N in [N : N in IdealsUpTo(2000, F) | not IsZero(N)] do
  H := HeckeCharacterGroup(N, [1,2,3]);
  D := DirichletGroup(N);
  psi_ct := 0;
  chi_ct := 0;
  for psi in Elements(D) do
    flag := true;
    for eps in epses do
      if psi(eps) ne &*[Sign(Evaluate(eps, v)) : v in vs] then
        flag := false;
      end if;
    end for;
    if flag then
      psi_ct +:= 1;
    end if;
  end for;
  for chi in Elements(H) do
    if IsCompatibleWeight(chi, k) then

      chi_ct +:= 1;
      print "---------------------------------";
      print "N", IdealOneLine(N);
      print "chi", chi, "order", Order(chi);
      print HeckeCharLabel(chi);
    end if;
  end for;
  assert chi_ct eq psi_ct;
end for;
