load "config.m";
import "ModFrmHil/definite.m" : BasisMatrixDefinite;

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


for N in IdealsUpTo(400, F) do
  pps := prime_exactly_div_level(N);
  H := HeckeCharacterGroup(N, [1,2]);
  chis := [chi : chi in Elements(H) | Order(chi) eq 1 and IsCompatibleWeight(chi, k)];
  for chi in chis do
    // check that there's a pp in pps such that 
    // the conductor of chi divides N / pp
    some_pp_works := (#pps gt 0) and &or[not (Conductor(chi) subset pp) : pp in pps]; 
    flag := false;
    if some_pp_works then
      print Norm(N), IdealOneLine(N);
      // find one that works
      for pp in pps do
        if not (Conductor(chi) subset pp) then
          flag := true;
          qq := pp;
          break;
        end if;
      end for;
      assert flag;
        
      // compute dimension using indefinite method
      B_indef := QuaternionAlgebra(qq, [InfinitePlaces(F)[1]]);
      O_indef := MaximalOrder(B_indef);
      Gamma := FuchsianGroup(O_indef);
      N_indef := N / qq;

      chi_indef := Restrict(chi, N / qq, [1,2]);
      dim_indef := ExactQuotient(Nrows(HeckeMatrix2(Gamma, N_indef, "Infinity", k, chi_indef)), 2);

      // compute dimension using definite method
      B_def := QuaternionAlgebra(1*ZF, InfinitePlaces(F) : Optimized);
      O_def := MaximalOrder(B_def);
      M := HilbertCuspForms(F, N, chi, k);
      dim_def := Nrows(BasisMatrixDefinite(M));

      print "dim_def:", dim_def, "dim indef:", dim_indef;
      if dim_def ne dim_indef then
        print "something's wrong!!";
      end if;
    end if;
  end for;
end for;
