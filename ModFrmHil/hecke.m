import !"Geometry/ModFrmHil/hecke.m" : make_ideal;

// For Fuchsian group Gamma, return a basis matrix for the plus subspace 
// of HeckeMatrix(Gamma, "Infinity") 

function basis_of_plus_subspace(M) 
  Gamma := FuchsianGroup(QuaternionOrder(M));
  N := Level(M) / make_ideal(Discriminant(QuaternionOrder(M)));
  T := HeckeMatrix2(Gamma, N, "Infinity");
  if T cmpeq [] then 
    T := Matrix(Integers(), 0, 0, []); 
  end if;
  T := ChangeRing(T, Integers());
  plus_basis := KernelMatrix(T-1);
  plus_basis := ChangeRing(plus_basis, Rationals()); 
  minus_basis := KernelMatrix(T+1);
  minus_basis := ChangeRing(minus_basis, Rationals()); 
  assert Nrows(plus_basis) + Nrows(minus_basis) eq Nrows(T);

  return plus_basis, minus_basis;
end function;

