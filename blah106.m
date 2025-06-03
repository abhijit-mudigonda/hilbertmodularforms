load "config.m";

SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

auts := Automorphisms(F);
aut := auts[2];

f := function(I)
  assert I subset 1*ZF;
  _, x := IsPrincipal(I);
  if ideal<ZF | aut(x)> eq I then
    return true;
  else
    return false;
  end if;
end function;

check_chi := function(N, chi)
  flag := true;
  for pp in PrimesUpTo(Norm(N), F : coprime_to:=Norm(N)*ZF) do
    // pp split in F/Q
    if IsPrime(Norm(pp)) and Norm(pp) ne 7 then
      if chi(Norm(pp)*ZF)^4 ne 1 then
        flag := false;
        break;
      end if;
    else
      // pp inert in F/Q
      if chi(pp)^2 ne 1 then
        flag := false;
        break;
      end if;
    end if;
  end for;
  if Order(chi) eq 2 then
    assert flag;
  end if;
  return flag;
end function;

k := [1,1,2];
ideals := [I : I in IdealsUpTo(10000, F) | not IsZero(I)];
for N in ideals do
  // printIdealOneLine(N);
  for chi in Elements(HeckeCharacterGroup(N, [1,2,3])) do
    if check_chi(N, chi) and IsCompatibleWeight(chi, k) then
      // print "-------------------------- Found one!-------------------------";
      print IdealOneLine(N), chi;
      // print "--------------------------------------------------------------";
      // print [Norm(tup[1]) : tup in Factorization(N)];
      if Order(chi) ne 2 then
        print "order isn't 2, instead it's", Order(chi);
      end if;
      if {7, 8, 27} subset SequenceToSet([Norm(tup[1]) : tup in Factorization(N)]) then
        print "ramified at 2, 3, and 7";
      elif {7, 8} subset SequenceToSet([Norm(tup[1]) : tup in Factorization(N)]) then
        print "ramified at 2 and 7";
      elif {8, 27} subset SequenceToSet([Norm(tup[1]) : tup in Factorization(N)]) then
        print "ramified at 2 and 3";
      elif {7, 27} subset SequenceToSet([Norm(tup[1]) : tup in Factorization(N)]) then
        print "ramified at 3 and 7";
      end if;
    end if;
  end for;
end for;
