import !"Geometry/ModFrmHil/hecke.m" : make_ideal;
import !"Geometry/ModFrmHil/indefinite.m" : ElementOfNormMinusOne;
import "weight_rep.m" : is_paritious;
import "hecke_field.m" : hecke_matrix_field;
import "diamond.m" : GetHeckeMatrix;

// For Fuchsian group Gamma, return a basis matrix for the plus subspace 
// of HeckeMatrix(Gamma, "Infinity") 

function basis_of_plus_subspace(M) 
  Gamma := FuchsianGroup(QuaternionOrder(M));
  N := Level(M) / make_ideal(Discriminant(QuaternionOrder(M)));
  assert N eq Level(M);
  T := GetHeckeMatrix(M, "Infinity" : SaveAndLoad:=true);
  // T := HeckeMatrix2(Gamma, N, "Infinity", Weight(M), DirichletCharacter(M));
  if T cmpeq [] then 
    T := Matrix(Integers(), 0, 0, []); 
  end if;
  if is_paritious(Weight(M)) then
    plus_basis := KernelMatrix(T-1);
    minus_basis := KernelMatrix(T+1);
  else
    // TODO abhijitm I don't think it matters which one is + and which is -,
    // since the other Hecke operators should act the same on both?
    T := ChangeRing(T, hecke_matrix_field(M));
    assert #Eigenvalues(T) eq 2;
    eig := Rep(Eigenvalues(T))[1];

    O := QuaternionOrder(M);
    F := BaseField(M);

    eps := Norm(ElementOfNormMinusOne(O));
    auts := AutsOfKReppingEmbeddingsOfF(F, F);
    k := Weight(M);
    k0 := Max(Weight(M));
    z := &*[auts[i](eps)^(k0 - k[i]) : i in [1 .. #auts]];

    assert F!(eig^2) eq z^-1;
    plus_basis := KernelMatrix(T - eig);
    minus_basis := KernelMatrix(T + eig);
  end if;
  assert Nrows(plus_basis) + Nrows(minus_basis) eq Nrows(T);
  assert Nrows(plus_basis) eq Nrows(minus_basis);

  plus_basis := ChangeRing(plus_basis, FieldOfFractions(BaseRing(plus_basis)));
  minus_basis := ChangeRing(minus_basis, FieldOfFractions(BaseRing(minus_basis)));
  
  return plus_basis, minus_basis;
end function;

