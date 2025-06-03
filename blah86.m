// code for Deding Yang

load "config.m";

F := QuadraticField(5);
ZF := Integers(F);

N := Factorization(41*ZF)[1][1];
H := HeckeCharacterGroup(N, [1,2]);
chi := H.1;

prec := 400;
M := GradedRingOfHMFs(F, prec);

for k in [[1, 9], [1, 7], [1, 5]] do
  Mk := HMFSpace(M, N, k, chi);
  Dk := DihedralBasis(Mk);

  assert #Dk eq 1;
  f := Dk[1];
  f := f / Coefficient(f, 1*ZF);

  _, mp := quo<ZF | 7*ZF>;
  for pp in PrimesUpTo(prec, F) do
    printf "%8o\t%25o\t%8o\n", LMFDBLabel(pp), Coefficient(f, pp), mp(Coefficient(f, pp));
  end for;
end for;
