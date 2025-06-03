load "config.m";


function prime_exactly_dividing_I(I)
  facts := Factorization(I);
  // check facts are sorted
  norms := [Norm(t[1]) : t in facts]; 
  assert norms eq Sort(norms);
  if exists(P){t[1] : t in facts | t[2] eq 1} then
    return true, P;
  else
    return false, _;
  end if;
end function;

// SetVerbose("ModFrmHil", 3);
// SetVerbose("HilbertModularForms", 3);

F := QuadraticField(5);
ZF := Integers(F);

k := [2,3];

N := ideal<ZF | 58, 4*ZF.2 + 34>;
print IdealOneLine(N);
b, qq := prime_exactly_dividing_I(N);
if b then
  primes := [pp : pp in PrimesUpTo(100, F) | IsCoprime(pp, N * qq)];
  B := QuaternionAlgebra(qq, [InfinitePlaces(F)[1]]);
  O := MaximalOrder(B);
  Gamma := FuchsianGroup(O);
  Ks := QuadraticExtensionsWithConductor(N, [1,2]);
  Ks := [K : K in Ks | IsTotallyComplex(K)];
  
  inert_primes_by_K := [];
  for K in Ks do
    ZK := Integers(K);
    inert_primes := {pp : pp in primes | IsPrime(ZK!!pp)};
    print "inert_primes size", #inert_primes;
    Append(~inert_primes_by_K, inert_primes);
  end for;

  H := HeckeCharacterGroup(N, [1,2]);
  N_JL := N / qq;
  for chi in Elements(H) do
    if (N_JL subset Conductor(chi)) and IsCompatibleWeight(chi, k) then
      chi_res := Restrict(chi, N_JL, [1,2]);
      hecke_matrices := AssociativeArray();
      for pp in primes do
        hecke_matrices[pp] := HeckeMatrix2(Gamma, N_JL, pp, k, chi_res);
      end for;
      if Nrows(hecke_matrices[primes[1]]) gt 0 then
        print "There are some [2,3] forms here!", IdealOneLine(N), chi, HeckeCharLabel(chi);
      end if;
      zero_eigs := {pp : pp in primes | Dimension(Kernel(hecke_matrices[pp])) gt 0};
      for i->inert_primes in inert_primes_by_K do
        if inert_primes subset zero_eigs then
          print "There might be a dihedral form here!", Ks[i], IdealOneLine(N), chi, HeckeCharLabel(chi);
        end if;
      end for;
    end if;
  end for;
end if;
