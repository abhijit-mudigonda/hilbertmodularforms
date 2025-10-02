load "config.m";

F := QuadraticField(5);
ZF := Integers(F);
M := GradedRingOfHMFs(F, 444);
k := [2, 2];

/*
N := ideal<ZF | 4*ZF.2 - 2>;
H := HeckeCharacterGroup(N, [1, 2]);
chi := H.1;
assert HeckeCharLabel(chi) eq "-5.0.1_20.1_2u1u1.2u";

Mk := HMFSpace(M, N, k, chi);
Sk := HeckeStabilityCuspBasis(Mk : prove:=false, stable_only:=true);
assert #Sk eq 2;
*/

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


for N in IdealsUpTo(100, F) do
  pps := prime_exactly_div_level(N);
  print Norm(N), IdealOneLine(N);
  H := HeckeCharacterGroup(N, [1,2]);
  chis := [chi : chi in Elements(H) | Order(chi) gt 2 and IsCompatibleWeight(chi, k)];
  for chi in chis do
    // check that there's a pp in pps such that 
    // the conductor of chi divides N / pp
    
    some_pp_works := (#pps gt 0) and &or[not (Conductor(chi) subset pp) : pp in pps]; 
    if some_pp_works then
      Mk := HMFSpace(M, N, k, chi);
      Sk := HeckeStabilityCuspBasis(Mk : prove:=false, stable_only:=true);
      if #Sk gt 0 then
        print "found one!", chi, HeckeCharLabel(chi), #Sk;
      end if;
  end for;
  end if;
end for;
