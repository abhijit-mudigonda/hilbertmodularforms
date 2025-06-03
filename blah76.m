load "config.m";

SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);
epses_F := UnitsGenerators(F);
vs := InfinitePlaces(F);

N := 2^3 * Factorization(7*ZF)[1][1];
H := HeckeCharacterGroup(N, [1,2,3]);
for chi in Elements(H) do
  if IsCompatibleWeight(chi, [1,1,2]) then
    print chi, Order(chi), IdealOneLine(Conductor(chi));
  end if;
end for;

K_abs := CyclotomicField(7);
IsSubfield(F, K_abs);
K := RelativeField(F, K_abs);
ZK := Integers(K);
ZK_abs := Integers(K_abs);
N_psi_lcm := ZK!!(N / Discriminant(ZK));
for N_psi_rel in Divisors(N_psi_lcm) do
  if N eq Norm(N_psi_rel) * Discriminant(ZK) then
    N_psi := ZK_abs!!N_psi_rel;
    D := DirichletGroup(N_psi);
    epses := UnitsGenerators(K_abs : exclude_torsion:=false);
    tors := epses[1];

    U, mU := UnitGroup(K_abs);
    assert Order(U.1) ne 0;
    for psi in Elements(D) do
      flag := true;
      b, order := IsRootOfUnity(psi(tors));
      assert b; // b should be a root of unity
      if order ne Order(U.1) then
        flag := false;
        continue;
      end if;
      for eps in epses_F do
        if psi(ZK_abs!eps) ne Sign(Evaluate(eps, vs[3])) then
          flag := false; 
        end if;
      end for;
      if flag then
        print "Found one", psi, IdealOneLine(Conductor(psi)), Norm(Conductor(psi));
      end if;
    end for;
  end if;
end for;

