load "config.m";

SetVerbose("ModFrmHil", 3);
SetVerbose("HilbertModularForms", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

function test_N_chi(N, chi, chi_label : k:=[1,1,2], prec:=200, chi_eis:=0)
  assert HeckeCharLabel(chi) eq chi_label;
  M := GradedRingOfHMFs(F, prec);
  Mk := HMFSpace(M, N, k, chi);
  Sk := HeckeStabilityCuspBasis(Mk : prove:=false, stable_only:=true, chi_eis:=chi_eis);

  if #Sk gt 0 then
    print "********* Found one!!!";
  end if;
  return Sk;
end function;


N := ideal<ZF | [182, 0, 0], [-80, -2, -2]>;
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.1;
chi_eis := H.1*H.2;
chi_label := "1.-2.-1.1_728.3_2u1.0u1.2.3u";
print "chi", Factorization(Conductor(chi));
print "N", Factorization(N);

print "chi 223", Factorization(Conductor(H.2));

test_N_chi(N, chi, chi_label : k:=[1,1,2], prec:=600, chi_eis:=chi_eis);
