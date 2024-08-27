import !"Geometry/ModFrmHil/hecke.m" : make_ideal;

// For Fuchsian group Gamma, return a basis matrix for the plus subspace 
// of HeckeMatrix(Gamma, "Infinity") 

function basis_of_plus_subspace(M) 
  Gamma := FuchsianGroup(QuaternionOrder(M));
  N := Level(M) / make_ideal(Discriminant(QuaternionOrder(M)));
  T := HeckeMatrix2(Gamma, N, "Infinity", Weight(M), DirichletCharacter(M));
  if T cmpeq [] then 
    T := Matrix(Integers(), 0, 0, []); 
  end if;
  plus_basis := KernelMatrix(T-1);
  minus_basis := KernelMatrix(T+1);
  assert Nrows(plus_basis) + Nrows(minus_basis) eq Nrows(T);

  plus_basis := ChangeRing(plus_basis, FieldOfFractions(BaseRing(plus_basis)));
  minus_basis := ChangeRing(minus_basis, FieldOfFractions(BaseRing(minus_basis)));
  
  return plus_basis, minus_basis;
end function;

