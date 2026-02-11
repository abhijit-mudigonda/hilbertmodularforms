Ds := [D : D in [1 .. 100] | IsSquarefree(D)];
Ds := [D : D in Ds | NarrowClassNumber(QuadraticField(D)) eq 1];
Ds := Ds[2 .. #Ds];

print Ds;

for D in Ds do
  F := QuadraticField(D);
  B := QuaternionAlgebra(1*Integers(F), InfinitePlaces(F));
  O := MaximalOrder(B);
  print D, #NarrowClassGroup(F), #RightIdealClasses(O);
end for;

