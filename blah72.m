load "config.m";

SetVerbose("ModFrmHil", 3);
SetVerbose("HilbertModularForms", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

function test_N_chi(N, chi : chi_label:=0, k:=[1,1,2], prec:=200)
  if chi_label cmpne 0 then
    assert HeckeCharLabel(chi) eq chi_label;
  end if;
  M := GradedRingOfHMFs(F, prec);
  Mk := HMFSpace(M, N, k, chi);
  Sk := HeckeStabilityCuspBasis(Mk : prove:=false, stable_only:=true);

  if #Sk gt 0 then
    print "********* Found one!!!";
  end if;
  return Sk;
end function;

N := 7*ZF;
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.1;
test_N_chi(N, chi : prec:=100);
