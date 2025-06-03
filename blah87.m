load "config.m";

SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

k := [1,1,2];
pp := Factorization(7*ZF)[1][1];
for t in [1 .. 22] do
  for s in [0 .. t] do
    for a in [0 .. s] do
      N := pp^a * 2^(s - a) * 3^(t-s);
      H := HeckeCharacterGroup(N, [1,2,3]);
      for chi in Elements(H) do
        if Order(chi) eq 2 or Order(chi) eq 6 then
          if IsCompatibleWeight(chi, k) then
            print a, s-a, t-s, Norm(N);
            break chi;
          end if;
        end if;
      end for;
    end for;
  end for;
end for;
