// find (N, chi) tuples with [2,2] forms and nebentypus order > 2

load "config.m";

F := QuadraticField(5);
ZF := Integers(F);
M := GradedRingOfHMFs(F, 300);
k := [2, 2];

function prime_exactly_div_level(I)
  facts := Factorization(I);
  norms := [Norm(t[1]) : t in facts]; 
  assert norms eq Sort(norms);
  pps_exactly_div := [];
  for tup in facts do
    if tup[2] eq 1 then
      Append(~pps_exactly_div, tup[1]);
    end if;
  end for;
  return pps_exactly_div;
end function;


for N in IdealsUpTo(200, F) do
  print Norm(N), IdealOneLine(N);
  H := HeckeCharacterGroup(N, [1,2]);
  chis := [chi : chi in Elements(H) | IsCompatibleWeight(chi, k) and Order(chi) gt 2];
  for chi in chis do
    Mk := HMFSpace(M, N, k, chi);
    if #HeckeStabilityCuspBasis(Mk : prove:=false, stable_only:=true) gt 0 then
      print "hi!!!", chi, Order(chi), HeckeCharLabel(chi);
    end if;
  end for;
end for;
