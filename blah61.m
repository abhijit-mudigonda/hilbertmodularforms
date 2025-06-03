load "config.m";

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);
epses_F := UnitsGenerators(F : exclude_torsion:=false);
vs := InfinitePlaces(F);

ideals := [N : N in IdealsUpTo(1000, F) | not IsZero(N)];
pp := ideal<ZF | [-7, 5, -2] >;
assert Norm(pp) eq 379;
for I in ideals do
  N := pp * I;
  print "N", IdealOneLine(N);
  Ks := QuadraticExtensionsWithConductor(N, [1,2,3]);
  Ks := [K : K in Ks | IsTotallyComplex(K)];
  for K in Ks do
    ZK := Integers(K);
    K_abs := AbsoluteField(K);
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
            print "Found one", psi, K, Norm(N);
          end if;
        end for;
      end if;
    end for;
  end for;
end for;
