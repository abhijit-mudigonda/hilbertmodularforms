load "config.m";

SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

prime_norms := [];

for p in PrimesUpTo(500) do
  if IsPrime(p*ZF) then
    Append(~prime_norms, p^3);
  else
    Append(~prime_norms, p);
  end if;
end for;

for N in IdealsUpTo(500, F) do
  if not (IsOne(N) or IsZero(N)) then
    D := DirichletGroup(N, [1,2,3]);
    for psi in Elements(D) do
      if &and[IsOne(psi(p)) : p in prime_norms] then
        print "found one!", psi, IdealOneLine(N);
      end if;
    end for;
  end if;
end for;
