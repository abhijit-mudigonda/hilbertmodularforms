load "config.m";

F := QuadraticField(5);
ZF := Integers(F);
M := GradedRingOfHMFs(F, 300);
N := 4*ZF;

ideals := [I : I in IdealsUpTo(100, F) | not IsZero(I)];

f := func<I | IdealToRep(M, I)>;
g := func<I | &*[IdealToRep(M, pair[1]^pair[2]) : pair in Factorization(I)]>;

/*
for I in ideals do
  facs := Factorization(I);
  nu := IdealToRep(M, I);
  mu := F!1;
  for pair in facs do
    mu *:= IdealToRep(M, pair[1]^pair[2]);
  end for;
  print "Norm(nu)", Norm(nu);
  print "Norm(mu)", Norm(mu);
end for;
*/
