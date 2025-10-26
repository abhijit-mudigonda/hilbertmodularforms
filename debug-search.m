load "config.m";
import "ModFrmHil/definite.m" : BasisMatrixDefinite, HeckeOperatorDefiniteBig;
import !"Geometry/ModFrmHil/hecke.m" : restriction;

verbose := false;

MAX_N_NORM := 300;

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

procedure test(F, k, neb_order : max_N_norm:=MAX_N_NORM)
  // F::FldAlg - A quadratic field
  // k::SeqEnum[RngIntElt] - A weight
  // neb_order::RngIntElt - Character orders to test

  assert Degree(F) eq 2;
  ZF := Integers(F);
  M := GradedRingOfHMFs(F, 100);

  indef_quat_data := AssociativeArray();

  for N in IdealsUpTo(max_N_norm, F) do
    pps := prime_exactly_div_level(N);
    H := HeckeCharacterGroup(N, [1,2]);
    chis := [chi : chi in Elements(H) | Order(chi) eq neb_order and IsCompatibleWeight(chi, k)];
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
        if not IsDefined(indef_quat_data, qq) then
          B_indef := QuaternionAlgebra(qq, [InfinitePlaces(F)[1]]);
          O_indef := MaximalOrder(B_indef);
          Gamma := FuchsianGroup(O_indef);
          indef_quat_data[qq] := <B_indef, O_indef, Gamma>;
        else
          B_indef, O_indef, Gamma := Explode(indef_quat_data[qq]);
        end if;

        N_indef := N / qq;

        _, K, _ := Splittings(B_indef);
        L := Compositum(K, CyclotomicField(Order(chi)));

        chi_indef := Restrict(chi, N / qq, [1,2]);
        primes := PrimesUpTo(100, F : coprime_to:=qq)[1 .. 3];
        dim_indef := ExactQuotient(Nrows(HeckeMatrix2(Gamma, N_indef, primes[1], k, chi_indef : HeckeMatrixField:=L)), 2);


        // compute dimension using definite method
        M := HilbertCuspForms(F, N, chi, k);
        // BasisMatrixDefinite by default will exclude Eisenstein series, so 
        // this is fine even in the trivial nebentypus case.
        dim_def := Nrows(BasisMatrixDefinite(M));

        M_qq_old := HilbertCuspForms(F, N / qq, chi_indef, k);
        dim_qq_old_def := Nrows(BasisMatrixDefinite(M_qq_old));

        // the indefinite method will see all the forms which are new at qq,
        // so we get the forms which are old at qq by looking at level N / qq.
        assert dim_def eq (dim_indef + 2 * dim_qq_old_def);

        /* 
        // TODO abhijitm, this is a further test to make sure that the big Hecke matrices are
        // in agreement. It's failing rn for some coefficient field reason.
        // primes should also probably include a prime dividing N.
        for pp in primes do
          def_hecke_mtrx := HeckeOperatorDefiniteBig(M, pp);
          def_qq_old_hecke_mtrx := HeckeOperatorDefiniteBig(M_qq_old, pp);
          def_hecke_mtrx := restriction(M, def_hecke_mtrx);
          def_qq_old_hecke_mtrx := restriction(M_qq_old, def_qq_old_hecke_mtrx);
          indef_hecke_mtrx := HeckeMatrix2(Gamma, N_indef, pp, k, chi_indef : HeckeMatrixField:=L);

          char_poly_full := CharacteristicPolynomial(def_hecke_mtrx);
          char_poly_qq_old := CharacteristicPolynomial(def_qq_old_hecke_mtrx);
          char_poly_qq_new := SquareRoot(CharacteristicPolynomial(indef_hecke_mtrx));
          assert char_poly_full eq (char_poly_qq_old)^2 * char_poly_qq_new;
        end for;
        */
      end if;
    end for;
  end for;
end procedure;

// TODO abhijitm test a field with narrow class number bigger than 1
Fs := [QuadraticField(5), QuadraticField(2)];

ks := [[2, 2], [3, 3], [2, 4], [2, 3]];

// TODO abhijitm get things to work with nebentypus order bigger than 2
neb_orders := [1, 2];

test(QuadraticField(5), [3, 3], 2);

/*
for F in Fs do
  for k in ks do
    for neb_order in neb_orders do
      test(F, k, neb_order);
    end for;
  end for;
end for;
*/
