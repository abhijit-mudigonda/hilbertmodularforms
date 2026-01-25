//freeze;

/**********************************************************************/
/*                                                                    */
/*  Non paritious Hilbert modular forms over Q(sqrt(2))               */
/*                                                                    */
/*           Created in November 2016                                 */
/*                                                                    */
/*                       By                                           */
/*                                                                    */
/*  Lassina Dembele, David Loeffler and Ariel Pacetti                 */
/*                                                                    */
/*  Last modified November 2016, Lassina Dembele                      */
/*                                                                    */
/**********************************************************************/


import "/Applications/Magma/package/Geometry/ModFrmHil/precompute.m": get_tps;


/**********************************************************************
  HilModFrm attributes
  Let M be a HilModFrm over the FldNum F
***********************************************************************/

declare type HilModFrm [HilModFrmElt];

declare attributes HilModFrm :

  BaseField,                 // F = BaseField(M), required to be an absolute field (for efficiency)
  Level,                     // Level(M) = an ideal in Integers(F)
  DirichletCharacter,        // always assigned: either the integer 1, or a GrpDrchNFElt with modulus Level(M)

  Weight,                    // Weight(M) = a sequence of integers, corresponding to InfinitePlaces(F)
  CentralCharacter,          // Only used internally.

  Dimension,                 // Dimension of the space of Hilbert modular forms.

  is_cuspidal,               // Is assigned when the space is cuspidal.

  HilModSpace,               // The space of Hilbert modular forms as created by Magma.
                             // We need this in order to call certain internal routines.

  QCD,                       // The data necessary for working with Hamilton quaternion algebra over Q(sqrt(5)).

  /////  BASIS MATRICES  /////

  basis_matrix,              // The rows generate M inside the vector space V used to compute it 
  basis_matrix_inv,          // basis_matrix * basis_matrix_inv = identity

  /////  THE STUFF NEEDED FOR THE PROJECTIVE LINE  /////

  PLD,                       // The projective line data required to work with P^(O/N).

  /////  THE WEIGHT REPRESENTATION  /////

  weight_rep,                // The weight representation: B^* -> GL_2(K)^g -> GL(V_k).

  weight_dimension,          // The dimension of the weight space.

  weight_base_field;         // The base field of the weight representation,
                             // or Rationals() for parallel weight 2


/*************************************************************************/


/* Here we define several record formats. */

ProjectiveLineData := recformat< FD : SeqEnum,            // This is a fundamental domain for the action of the
                                                          // the unit group in some quaternion algebra on P^1(O/n).

                                 orbits : SeqEnum,        // The sequence containing the orbits the quaternion action
                                                          // on P^1(O/n).

                                 Stabs : SeqEnum,         // Generators of the group of stabilizers, for each element
                                                          // in FD.

                                 StabOrders : SeqEnum,    // Order of the group of stabilizers, for each element
                                                          // in FD.

                                 Lookuptable : Assoc,     // The lookuptable for the elements in P1List.

                                 P1List : SetIndx,        // Elements of the projective line P^1(O/n).

                                 P1Rep : UserProgram,     // Function sending [u,v] to an element of P1List.

                                 splitting_map : Map,     // The splitting map mod n for the quaternion order.

                                 Level : RngOrdIdl        // The level.
                               >;


QuaternionicData := recformat< QuaternionOrder : AlgAssVOrd, // Quaternion order used for the computation.

                               max_order_units : SeqEnum,    // The unit group in the quaternion order.

                               unit_map : Map                // The map from the unit group of the quaternion
                             >;                              // as an abstract group to the quaternions.


/*************************************************************************/


/*************************************************************************
  Printing for HilModFrm
**************************************************************************/


function ideal_print_on_one_line(I, level)
  if Type(I) in {RngIntElt,FldRatElt} then
    printf "%o", I; return 0;
  elif Type(I) eq RngInt then 
    printf "%o", Minimum(I); return 0; 
  end if;
  gens := Generators(I);
  printf "Ideal of norm %o generated by ", Norm(I);
  for g := 1 to #gens-1 do printf "%o, ", gens[g]; end for;
  printf "%o", gens[#gens];
  return 0;
end function;

intrinsic Print(x::HilModFrm, level::MonStgElt)
{}
  Print(x, level);
end intrinsic;

intrinsic Print(x::HilModFrm, level::MonStgElt)
{}
  F := BaseField(x);

  assert assigned x`is_cuspidal;

  if x`is_cuspidal and IsNew(x) then
    printf "Cuspidal new subspace of Hilbert modular forms";
  elif x`is_cuspidal and Norm(Level(x)) ne 1 then
    printf "Cuspidal subspace of Hilbert modular forms";
  elif x`is_cuspidal then
    printf "Cuspidal subspace of Hilbert modular forms";
  else 
    printf "Space of Hilbert modular forms";
  end if;

  if level ne "Minimal" then 
    printf " over";
    printf "\n    %o", BaseField(x);
    printf "\n    Level: "; 
    _ := ideal_print_on_one_line(Level(x), level);
    if x`DirichletCharacter cmpne 1 then
      printf "\n    Hecke character: %o", x`DirichletCharacter;
    end if;
    printf "\n    Weight: %o", Weight(x);
    if assigned x`Dimension then 
      printf "\n    Dimension %o", x`Dimension; 
    end if;
  end if;
  if level eq "Maximal" and assigned x`QuaternionOrder then 
    OM := x`QuaternionOrder;
    printf "\ncomputed as automorphic forms on %o order in definite quaternion algebra of discriminant ",
           IsMaximal(OM) select "maximal" else "Eichler";
    _ := ideal_print_on_one_line( Discriminant(Algebra(OM)), level);
  end if;
end intrinsic;


intrinsic IsCoercible(x::HilModFrm, y::.) -> BoolElt, .
{}
  if Type(y) eq HilModFrmElt and Parent(y) cmpeq x then
    return true, y;
  end if;
  return false;
end intrinsic;

intrinsic 'in'(x::., y::HilModFrm) -> BoolElt
{}
  if Type(x) ne HilModFrmElt then
    return false, "The first argument should be a HilModFrmElt";
  else
    return Parent(x) eq y;
  end if;
end intrinsic;




intrinsic BaseField(M::HilModFrm) -> Fld
  {The base field of the space of Hilbert modular forms}
  return M`BaseField;
end intrinsic;





/**********************************************************************
  Attributes of an HilModFrmElt...
  Let f be a HilModFrmElt, contained in the HilModFrm space M

***********************************************************************/

declare type HilModFrmElt;

declare attributes HilModFrmElt :
  
  Parent,                       // = M

  BaseField,                    // the field over which the element is defined,

  coords_wrt_parent,            // Vector representing element wrt the rational basis of M
  coords_wrt_ambient,           //  ..... of M`ambient
  coords_wrt_top,               //  ..... of top space in chain M`ambient`ambient
                                //        [mercifully, not using the top one yet]
  coords_raw,                   // = coords_wrt_parent * M`basis_matrix

  hecke_eigenvalues;            // array indexed by primes


/****************   Hacking for HilModFrmElt   *****************/

intrinsic Parent(x::HilModFrmElt) -> HilModFrm
{}
  return x`Parent;
end intrinsic;


intrinsic Print(x::HilModFrmElt, level::MonStgElt)
{}
  if level eq "Maximal" then 
    printf "Element of %o\n", Parent(x); 
    printf "defined over %o", BaseField(x); 
  else
    printf "Element of "; 
    Print(Parent(x), "Minimal");
  end if; 
end intrinsic;


intrinsic IsCoercible(x::HilModFrmElt, y::.) -> BoolElt, .
{}
  return false;
end intrinsic;

intrinsic 'in'(x::., y::HilModFrmElt) -> BoolElt
{}
  return false;
end intrinsic;


intrinsic BaseField(f::HilModFrmElt) -> Fld
{The field over which the Hilbert modular form f is defined}
  return f`BaseField;
end intrinsic;







/**********************************************************************************/



/* Here we define the routines needed for the weight representation .*/


intrinsic HidaNormalization(k::SeqEnum[RngIntElt]) -> RngIntElt, SeqEnum, SeqEnum
  { Given a sequence of integers k, view as a weight for Hilbert modular forms, this returns
    a positive integer k_0, and two sequences m and n, where m is one of half integers m such
    that k_i + 2*m_i is independent of i, and n_i = k_i - 2. }

  require forall(i){ i : i in [1..#k] | k[i] ge 2 } : "The weight must be cohomological!";

  k0 := Max(k);
  m := [ (k0 - k[i])/2 : i in [1..#k] ];
  n := [ k[i] - 2 : i in [1..#k] ];
  return k0, m, n;
end intrinsic;

      

intrinsic IsParitious(k::SeqEnum[RngIntElt]) -> BoolElt
  { Given a sequence of non negative integers, i.e. a weight for Hilbert modular forms,
    this returns true if k is paritious, and false otherwise. }

  _, m := HidaNormalization(k);
  if m subset Integers() then
    return true;
  else
    return false;
  end if;
end intrinsic;

      

// Tensor product is associative; for efficiency always do
// TensorProduct(small_mat, large_mat)

function arithmetic_weight_map(q, splittings, K1, m, n)
  d := #m;
  m1 := [ Integers() | s : s in m ];
  M1 := MatrixRing(K1, 1)!1;
  for l := d to 1 by -1 do
    if m1[l] eq 0 and n[l] eq 0 then
      // don't need to modify M
      continue;
    else
      matq := splittings[l](q);
      if n[l] eq 0 then
        M1 *:= Determinant(matq)^m1[l];
      else
        if n[l] eq 1 then
          Ml := matq;
        else
          Ml := SymmetricPower2(matq, n[l]);
        end if;
        if m[l] ne 0 then
          Ml *:= Determinant(matq)^m1[l];
        end if;
        if l eq d then
          M1 := Ml;
        else
          M1 := TensorProduct(Ml, M1);
        end if;
      end if;
    end if;
  end for;
  return M1;
end function;

      
//Here, we might need to take square roots as the sequence m is
// one of possibly half integers.

function non_arithmetic_weight_map(q, splittings, K1, n)
  d := #n;
      
  M1 := MatrixRing(K1, 1)!1;
  mats := [ splittings[l](q) : l in [1..d] ];
      
  //Evaluating the weight representation.

  for l := 1 to d do
    if n[l] eq 0 then
      // don't need to modify M
      continue;
    else
      if n[l] eq 1 then
        Ml := mats[l];
      else
        Ml := SymmetricPower2(mats[l], n[l]);
      end if;
      M1 := TensorProduct(Ml, M1);
    end if;
  end for;
  return M1;
end function;


      

function weight_representation(M)
// Given a space of Hilbert modular forms over Q(sqrt(5)). This returns the weight representation.
// This code is specific to F = Q(sqrt(5)) in the sense that it uses K = Q(zeta_5) as the Galois
// closure of F over Q where H is splits. The function returns a map H^* -> GL(2, K)^2 -> GL(V_k),
// where V_k the weight space.

  if assigned M`weight_rep then
    return M`weight_rep, M`weight_dimension, M`weight_base_field;
  else
    OH := M`QCD`QuaternionOrder;
    H<i1,j1,k1> := Algebra(OH);
    F := BaseField(H);
    O<w> := Integers(F);
    k := M`Weight;
    k0, m, n := HidaNormalization(k);

    if Seqset(k) eq {2} then // parallel weight 2
      I := IdentityMatrix(Rationals(), 1);
      Mat1 := Parent(I);
      M`weight_rep := map< H -> Mat1 | q :-> I >;
      M`weight_base_field := Rationals();
      M`weight_dimension := 1;
    else
      K<z> := CyclotomicField(8);
      OK := Integers(K);
      pol := DefiningPolynomial(F);
      rts := Roots(pol, K);
      embeddings_F_to_K := [ hom< F -> K | r[1]> : r in rts ];

      // Archimedian splitting
      HK, toHK := ChangeRing(H, K);
      BH := [ 1, (1 + i1)/w, (1 + j1)/w, (1 + i1)*(1 + j1)/2 ];
      OHK := MaximalOrder(Order([ toHK(x) : x in BH ]));
      _, M2K, toM2K := IsMatrixRing(HK : Isomorphism := true);
      bool, _, t := IsConjugate(OHK, Order([ a@@toM2K : a in Basis(MatrixRing(OK, 2)) ]) : FindElement := true);
      BHK := [ toM2K(t*toHK(x)*t^-1): x in BH ];
      standard_basis_image := [ toM2K(t*toHK(x)*t^-1) : x in [1, i1, j1, k1 ] ];
      
      // The weight_base_field should be an absolute field for efficiency in later calculations
      M`weight_base_field := K;
      splitting_seq := [];
      for m1 := 1 to Degree(F) do
        h := embeddings_F_to_K[m1];
        iK := standard_basis_image[2];
        jK := standard_basis_image[3];
        kK := standard_basis_image[4];
        assert kK eq iK*jK;
        Append(~splitting_seq,
           map< H -> MatrixRing(K, 2) | q :-> h(s[1]) + h(s[2])*iK + h(s[3])*jK + h(s[4])*kK where s := Eltseq(q) >);
      end for;
      M`weight_dimension := &*[ x + 1 : x in n ];
      VK := MatrixRing(K, M`weight_dimension);
      if IsArithmeticWeight(F, k) then
        M`weight_rep := map< H -> VK | q :-> arithmetic_weight_map(q, splitting_seq, K, m, n) >;
      else
        M`weight_rep := map< H -> VK | q :-> non_arithmetic_weight_map(q, splitting_seq, K, n) >;
      end if;
    end if;
  return M`weight_rep, M`weight_dimension, M`weight_base_field;
  end if;
end function;
      
      


/* Computing the space of quaternionic automorphic forms. */

function minimal_generators(G)
  gens0 := Generators(G);
  gens := [];
  S := sub<G | >;
  for g in gens0 do
    if g notin S then
      Append(~gens, g);
      S := sub<G | gens>;
      if S eq G then
        break g;
      end if;
    end if;
  end for;
return gens;
end function;

     
      
function InvariantSpace(stab, weight_map, k, eps)
//Given a stabilizer stab (as sequence of quaternions), and a weight representation, this returns
//the corresponding space of invariant space. This will be used in the computation of the space of
//Hilbert modular forms.

  H1 := Domain(weight_map);
  WM := Codomain(weight_map);
  K1 := BaseRing(WM);
  weight_dim := Degree(WM);
  V := VectorSpace(K1, weight_dim);
  V1 := V;
  U, Umap := UnitGroup(BaseField(H1));
  if IsParitious(k) then
    S := [ [* H1!u@Umap, u@Umap *] : u in Generators(U) ];
    S := S cat stab;
  else
    S := stab cat [ [* H1!(-1), -1 *] ];
  end if;
      
  for s in S do
    g0 := hom< V -> V | weight_map(s[1]) - eps(s[2])^-1*ScalarMatrix(K1, weight_dim, 1)>;
    V1 := V1 meet Kernel(g0);
    if Rank(V1) eq 0 then
      break;
    end if;
  end for;
      
  V1 := sub< V | V1 >;
  return BasisMatrix(V1);
end function;



function concatenate_matrices(A, B);
  //Given to matrices A and B, this returns the block matrix [A, 0, 0, B].

  if Nrows(A) eq 0 then
    return B;
  elif Nrows(B) eq 0 then
    return A;
  else
    Q := ZeroMatrix(BaseRing(A), Nrows(A) + Nrows(B), Ncols(A) + Ncols(B));
    Q := InsertBlock(Q, A, 1, 1);
    Q := InsertBlock(Q, B, Nrows(A) + 1, Ncols(A) + 1);
    return Q;
  end if;
end function;



function HilbertModularSpaceBasisMatrix(M)
// Given a space of quaternion Hilbert modular forms, this returns the data relative to the projective line for the
// computation of the Brandt module. The output consists of the projective line along its orbits of the action of
// the units in the order $OH$.

d := M`Level;
k := M`Weight;
units := M`QCD`max_order_units;
unit_map := M`QCD`unit_map;
split_map := M`PLD`splitting_map;
weight_map := M`weight_rep;
eps := M`DirichletCharacter;

OK := Order(d);
R := quo<OK|d>;
P1, P1rep := ProjectiveLine(R: Type := "Matrix");
WM := Codomain(weight_map);
K1 := BaseRing(WM);
basis_matrix := ZeroMatrix(K1, 0, 0);

Lookuptable := AssociativeArray(P1);
FD := [];
StabOrders := [Integers()| ];
Stabs := [ ];

U := Domain(unit_map);
UU := Universe(units);
assert Type(UU) eq AlgQuat;
assert units[1] in { UU | 1, -1 };
mats := [ split_map(q) : q in units ];
elts := Set(P1);
while not IsEmpty(elts) do
  e := Random(elts);
  orbit := [];
  stab := [];
  elt_datas := [];
    
  for j := 1 to #units do
    e0 := mats[j] * e;
    _, e1, s1 := P1rep(e0, false, true);
    if e1 notin orbit then
      Exclude(~elts, e1);
      Append(~orbit, e1);
      Append(~elt_datas, [* j, s1 *]);
    elif e1 eq e then
      Append(~stab, [* units[j], s1 *]);
    end if;
  end for;
      
  //Here we compute generators for the stabilizer.
  stab_U := [s[1] @@ unit_map : s in stab];
  stab_U_gens := minimal_generators(sub< U | stab_U >);
  stab_gens := [stab[i] : i in [1..#stab] | stab_U[i] in stab_U_gens];
      
  //Determining the base matrix.
  invts_space := InvariantSpace(stab, weight_map, k, eps);
  if Rank(invts_space) ge 1 then
    Append(~FD, e);
    Append(~StabOrders, 1+#stab);
    Append(~Stabs, stab_gens);
    basis_matrix := concatenate_matrices(basis_matrix, invts_space);
    for l in [1..#orbit] do
      Lookuptable[orbit[l]] := Insert(elt_datas[l], 1, #FD);
    end for;
  end if;
end while;

if Rank(basis_matrix) ge 1 then
  basis_matrix := EchelonForm(basis_matrix);
  basis_matrix_inv := Transpose(Solution(Transpose(basis_matrix), ScalarMatrix(K1, Rank(basis_matrix), 1)));
else
  basis_matrix := ZeroMatrix(K1, 0, 0);
  basis_matrix_inv := ZeroMatrix(K1, 0, 0);
end if;
       
//We now create the record for the projective data.
PLD := rec< ProjectiveLineData |
            Level := d,
            P1List := P1,
            P1Rep := P1rep,
            FD := FD,
            splitting_map := split_map,
            Lookuptable := Lookuptable,
            StabOrders := StabOrders,
            Stabs := Stabs
          >;
M`PLD := PLD;
M`basis_matrix := basis_matrix;
M`basis_matrix_inv := basis_matrix_inv;
return M`PLD, M`basis_matrix, M`basis_matrix_inv;
end function;



function DivisorHilbertModularSpace(M, d)
// Given a divisor $d$ of the level, this creates the space of automorphic form of level d,
// and weight k.

assert GCD(d, M`Level) eq d;

N := New(HilModFrm);
N`Level := d;
N`Weight := M`Weight;
N`QCD := M`QCD;
N`weight_rep := M`weight_rep;
N`DirichletCharacter := Restrict(M`DirichletCharacter, d);
if Norm(d) eq 1 then
  N`PLD, N`basis_matrix, N`basis_matrix_inv := HilbertModularSpaceBasisMatrix(N);
  return N;
end if;

k := M`Weight;
units := (M`PCD)`max_order_units;
unit_map := (M`PCD)`unit_map;
split_map := (M`PLD)`splitting_map;
eps := Restrict(M`DirichletCharacter, d);
weight_map := M`weight_rep;
WM := Codomain(weight_map);
K1 := BaseRing(WM);
OK := Order(d);
R := quo<OK|d>;
P1, P1rep := ProjectiveLine(R: Type := "Matrix");
Lookuptable := AssociativeArray(P1);
FD := (M`PLD)`FD;
d_FD := [];
d_StabOrders := [Integers()| ];
d_Stabs := [ ];
basis_matrix := ZeroMatrix(K1, 0, 0);

U := Domain(unit_map);
UU := Universe(units);
assert Type(UU) eq AlgQuat;
assert units[1] in { UU | 1, -1 };
mats := [ split_map(q) : q in units ];
elts := { x1 where _, x1 := P1rep(x, false, false) : x in FD };

while not IsEmpty(elts) do
  e := Random(elts);
  orbit := [];
  stab := [];
  elt_datas := [];
  for j := 1 to #units do
    e0 := mats[j] * e;
    _, e1, s1 := P1rep(e0, false, true);
    if e1 notin orbit then
      Exclude(~elts, e1);
      Append(~orbit, e1);
      Append(~elt_datas, [* j, s1 *]);
    elif e1 eq e then
      Append(~stab, [* units[j], s1 *]);
    end if;
  end for;

  //Here we compute generators for the stabilizer.
  stab_U := [s[1] @@ unit_map : s in stab];
  stab_U_gens := minimal_generators(sub< U | stab_U >);
  stab_gens := [stab[i] : i in [1..#stab] | stab_U[i] in stab_U_gens];

  //Determining the base matrix.
  invts_space := InvariantSpace(stab, weight_map, k, eps);
  if Rank(invts_space) ge 1 then
    Append(~d_FD, e);
    Append(~d_StabOrders, 1+#stab);
    Append(~d_Stabs, stab_gens);
    basis_matrix := concatenate_matrices(basis_matrix, invts_space);
    for l in [1..#orbit] do
      Lookuptable[orbit[l]] := Insert(elt_datas[l], 1, #d_FD);
    end for;
  end if;
end while;

if Rank(basis_matrix) ge 1 then
  basis_matrix := EchelonForm(basis_matrix);
  basis_matrix_inv := Transpose(Solution(Transpose(basis_matrix), ScalarMatrix(K1, Rank(basis_matrix), 1)));
else
  basis_matrix := ZeroMatrix(K1, 0, 0);
  basis_matrix_inv := ZeroMatrix(K1, 0, 0);
end if;
       
//We now create the record for the projective data.
d_PLD := rec< ProjectiveLineData |
              Level := d,
              P1List := P1,
              P1Rep := P1rep,
              FD := d_FD,
              Lookuptable := Lookuptable,
              StabOrders := d_StabOrders,
              Stabs := d_Stabs,
              splitting_map := split_map
            >;
N`PLD := d_PLD;
N`basis_matrix := basis_matrix;
N`basis_matrix_inv := basis_matrix_inv;
return N;
end function;




/* Computing the degeneracy maps. */

function degeneracy_map1(M1, M2)
//Given two spaces Hilbert forms M1, M2 such that Level(M1) divides Level(M2)
//this returns the degeneracy map corresponding to the usual pull back
//P^1(Level(M2)) -> P^1(Level(M1)).
        
if Dimension(M1) eq 0 then
  return ZeroMatrix(M1`weight_base_field, 0, 0);
end if;

eps := M1`DirichletCharacter;
R := Domain(eps);
K1 := M1`weight_base_field;
PLD1 := M1`PLD;
FD1 := PLD1`FD;
FD2 := (M2`PLD)`FD;
P1 := PLD1`P1List;
P1rep := PLD1`P1Rep;
split_map := PLD1`splitting_map;
Lookuptable1 := PLD1`Lookuptable;
weight_map := M1`weight_rep;
units := M1`QCD`max_order_units;
wt_dim := M1`weight_dimension;
WM := Codomain(weight_map);
tp_mat := RMatrixSpace(K1, Ncols(M1`basis_matrix), Ncols(M2`basis_matrix))!0;
 
for m := 1 to #FD2 do
  _, u0, v0 := P1rep(FD2[m], true, true);
  if IsDefined(Lookuptable1, u0) then
    elt_data := Lookuptable1[u0];
    if not IsTrivial(eps) then
      v1 := (R!v0)^-1*(R!elt_data[3]);
      chi_val := eps(v1);
    else
      chi_val := K1!1;
    end if;
    wm := weight_map(units[elt_data[2]])^-1;

    mat0 := chi_val*wm;
    mat1 := ExtractBlock(tp_mat, (elt_data[1] - 1)*wt_dim + 1, (m - 1)*wt_dim + 1, wt_dim, wt_dim);
    mat1 := mat1 + mat0;
    tp_mat := InsertBlock(tp_mat, mat1, (elt_data[1] - 1)*wt_dim + 1, (m - 1)*wt_dim + 1);
  end if;
end for;
brandtp := M1`basis_matrix*tp_mat*M2`basis_matrix_inv;
       
return brandtp;
end function;


function degeneracy_map2(M1, M2, p, tp)
//Given two spaces Hilbert forms M1, M2, an element v in M1 and an ideal p such that Level(M1) divides Level(M2)
//and p divides Level(M2)/Level(M1), this returns the image of v in M2 under the degeneracy map.

if Dimension(M1) eq 0 then
  return ZeroMatrix(M1`weight_base_field, 0, 0);
end if;

d := M1`Level;
eps := M1`DirichletCharacter;
OK := Order(d);
R := quo<OK|d>;
K1 := M1`weight_base_field;
PLD1 := M1`PLD;
FD1 := PLD1`FD;
FD2 := (M2`PLD)`FD;
P1 := PLD1`P1List;
P1rep := PLD1`P1Rep;
split_map := PLD1`splitting_map;
Lookuptable1 := PLD1`Lookuptable;
weight_map := M1`weight_rep;
units := M1`QCD`max_order_units;
wt_dim := M1`weight_dimension;
WM := Codomain(weight_map);
tp_mat := RMatrixSpace(K1, Ncols(M1`basis_matrix), Ncols(M2`basis_matrix))!0;

q_seq := [ Norm(q)^-1*q : q in tp ];
mats := [ split_map(q) : q in tp ];
nr_seq := [ Norm(q)^-1 : q in tp ];

for l in [1..#mats], m in [1..#FD2] do
  u1 := Eltseq(mats[l]*FD2[m]);
  vals := [Valuation(s, p): s in u1];
  if Min(vals) ge 1 then
    if Valuation(d, p) eq 0 then
      uv := Matrix(OK, 2, 1, [OK!R!(nr_seq[l]*u1[1]), OK!R!(nr_seq[l]*u1[2])]);
    elif Min(vals) eq 1 then
      uv := Matrix(OK, 2, 1, [OK!R!(nr_seq[l]*u1[1]), OK!R!(nr_seq[l]*u1[2])]);
    end if;
    t0, u0, v0 := P1rep(uv, true, true);
    if t0 and IsDefined(Lookuptable1, u0) then
      elt_data := Lookuptable1[u0];
      if not IsTrivial(eps) then
        v1 := (R!v0)^-1*(R!elt_data[3]);
        chi_val := eps(v1);
      else
        chi_val := K1!1;
      end if;
      wm := weight_map(units[elt_data[2]]^-1*q_seq[l]);
      mat0 := chi_val*wm;
      mat1 := ExtractBlock(tp_mat, (elt_data[1] - 1)*wt_dim + 1, (m - 1)*wt_dim + 1, wt_dim, wt_dim);
      mat1 := mat1 + mat0;
      tp_mat := InsertBlock(tp_mat, mat1, (elt_data[1] - 1)*wt_dim + 1, (m - 1)*wt_dim + 1);
    end if;
  end if;
end for;
brandtp := M1`basis_matrix*tp_mat*M2`basis_matrix_inv;
return brandtp;
end function;



/* Computing Hecke operators. */

function hecke_operator(M, tp)

if Dimension(M) eq 0 then
  return ZeroMatrix(M`weight_base_field, Dimension(M), Dimension(M));
end if;

eps := M`DirichletCharacter;
R := Domain(eps);
K1 := M`weight_base_field;
PLD := M`PLD;
FD := PLD`FD;
P1 := PLD`P1List;
P1rep := PLD`P1Rep;
split_map := PLD`splitting_map;
Lookuptable := PLD`Lookuptable;
weight_map := M`weight_rep;
units := M`QCD`max_order_units;
wt_dim := M`weight_dimension;
WM := Codomain(weight_map);
A1 := M`basis_matrix;
tp_mat := ScalarMatrix(K1, Ncols(A1), 0);

mats := [ split_map(q) : q in tp ];
for l := 1 to #mats do
  for m := 1 to #FD do
    t0, u0, v0 := P1rep(mats[l]*FD[m], true, true);
    if t0 then
      if IsDefined(Lookuptable, u0) then
        elt_data := Lookuptable[u0];
        if not IsTrivial(eps) then
          v1 := (R!v0)^-1*(R!elt_data[3]);
          chi_val := eps(v1);
        else
          chi_val := K1!1;
        end if;
        wm := weight_map(units[elt_data[2]]^-1*tp[l]);
        mat0 := chi_val*wm;
        mat1 := ExtractBlock(tp_mat, (elt_data[1] - 1)*wt_dim + 1, (m - 1)*wt_dim + 1, wt_dim, wt_dim);
        mat1 := mat1 + mat0;
        tp_mat := InsertBlock(tp_mat, mat1, (elt_data[1] - 1)*wt_dim + 1, (m-1)*wt_dim + 1);
      end if;
    end if;
  end for;
end for;
    
return M`basis_matrix*tp_mat*M`basis_matrix_inv;
end function;



      
/* Creating the space of quaternionic Hilbert modular forms. */
      
intrinsic HilbertModularForms(N::RngOrdIdl, k::SeqEnum, eps::GrpDrchNFElt) -> HilModFrm
  { Given an integral ideal N in F = Q(sqrt(5)), a weight k and a Hecke character, this
    creates the space of automorphic forms of weight k, level $N$ and Hecke character eps over F. }

  M := New(HilModFrm);
  M`Level := N;
  M`Weight := k;
  O<w> := Order(N);
  F := NumberField(O);
  M`BaseField := F;
  M`DirichletCharacter := eps;

  // Creating the quaternion algebra
  H<i1,j1,k1> := QuaternionAlgebra< F | -1, -1 >;
  BH := [ 1, (1 + i1)/w, (1 + j1)/w, (1 + i1)*(1 + j1)/2 ];
  OH := MaximalOrder(Order(BH));
  UH, toOH := UnitGroup(OH);
  units := [ H!toOH(u) : u in UH ];

  // Quaternionic data
  QCD := rec< QuaternionicData | QuaternionOrder := OH,
                                 max_order_units := units,
                                 unit_map := toOH
            >;
  M`QCD := QCD;
  M`HilModSpace := HilbertCuspForms(F, 1*O : QuaternionOrder := OH);
      
  // Creating the projective line data
  _, split_map := ResidueMatrixRing(OH, N);
  PLD := rec< ProjectiveLineData | splitting_map := split_map >;
  M`PLD := PLD;

  // Creating the weight representation
  _ := weight_representation(M);

  return M;
end intrinsic;
      
      

      
intrinsic HilbertModularForms(N::RngOrdIdl, k::SeqEnum) -> HilModFrm
  { Given an integral ideal N in F = Q(sqrt(5)) and a weight k, this creates the space of automorphic
    forms of weight k, level $N$ and trivial Hecke character over F. }

  G := DirichletGroup(N);
  eps0 := G!1;
  return HilbertModularForms(N, k, eps0);
end intrinsic;
      



intrinsic HilbertModularForms(N::RngOrdIdl, eps::GrpDrchNFElt) -> HilModFrm
  { Given an integral ideal N in F = Q(sqrt(5)), and a Hecke character eps, this creates the space of automorphic
    forms of weight 2, level $N$ and Hecke character eps over F. }
      
  O := Order(N);
  k0 := [ 2 : i in [1..Degree(O)] ];
  return HilbertModularForms(N, k0, eps);
end intrinsic;
      



intrinsic Dimension(M::HilModFrm) -> RngIntElt
  { Given a space of Hilbert modular forms, this returns its dimension. }

  if not assigned M`Dimension then
    // The projective line data
    k := M`Weight;
    F := BaseField(M);
    O := Integers(F);
    eps := M`DirichletCharacter;
    weight_map := M`weight_rep;
    K1 := M`weight_base_field;
      
    if eps(-1) eq  weight_map(-1)[1,1] then
      _ := HilbertModularSpaceBasisMatrix(M);
    else
      M`basis_matrix := ZeroMatrix(K1, 0, 0);
      M`basis_matrix_inv := ZeroMatrix(K1, 0, 0);
    end if;
  end if;

  M`Dimension := Rank(M`basis_matrix);
  return M`Dimension;
end intrinsic;

      
      

/* Computing the degeneracy maps and the new subspace. */
      
intrinsic DegeneracyMap(M1::HilModFrm, M2::HilModFrm, p::RngOrdIdl) -> Mtrx
  {Given two spaces quaternionic automorphic forms M1, M2, such that Level(M2)
   divides Level(M2), and a prime ideal p which divides Level(M2)/Level(M1),
   this returns the downward degeneracy map M1 -> M2. This only works when
   DegeneracyMapDomain is used to create M2 from M1. The output is a matrix
   which represents the degeneracy map in the bases of M1 and M2.}

  N1 := M1`Level;
  N2 := M2`Level;
  require GCD(N1, N2) eq N2 : "Level(M1) must divide Level(M2)!";
    
  O := Order(N1);
  d := O!!(N1/N2);
  cond := Conductor(M2`DirichletCharacter);
  require GCD(p, d) eq p : "p must divide Level(M1)/Level(M2)!";
  require GCD(cond, N2) eq cond : "The conductor of the character must divide the level of M2!";

  MM := M1`HilModSpace;
  O := Integers(BaseField(M1));
  OH := M1`QCD`QuaternionOrder;
  H := Algebra(OH);
  units := M1`QCD`max_order_units;
  _, u := IsNarrowlyPrincipal(p);
  tp := get_tps(MM, u*O);
  return degeneracy_map2(M1, M2, p, tp);
end intrinsic;

      
      
      
intrinsic DegeneracyMapDomain(M::HilModFrm, p::RngOrdIdl) -> HilModFrm
  {Given a space of Hilbert modular forms and a prime p that divides the level of M,
    this creates a space of Hilbert modular forms of level Level(M)/p, whose data is
    compatible with M for the computation of the degeneracy maps.}

  require IsPrime(p) and GCD(M`Level, p) eq p :
                               "The second argument must be a prime that divides the level of the space M!";
  O := Order(p);
  d := O!!(M`Level/p);
  cond := Conductor(M`DirichletCharacter);

  if GCD(d, cond) ne cond then
    N := New(HilModFrm);
    N`basis_matrix := ZeroMatrix(M`weight_base_field, 0, 0);
    N`basis_matrix_inv := ZeroMatrix(M`weight_base_field, 0, 0);
  else
    N := DivisorHilbertModularSpace(M, d);
  end if;

  return N;
end intrinsic;



intrinsic NewSubspace(M::HilModFrm) -> HilModFrm
  {Given a space of Hilbert modular forms, this computes the subpace of newforms.}

  if Norm(M`Level) eq 1 or Dimension(M) eq 0 then
    if not assigned M`is_new then
      M`is_new := true;
    end if;
    return M;
  end if;
       
  K1 := M`weight_base_field;
  levelf := Factorization(M`Level);
  V := VectorSpace(K1, Rank(M`basis_matrix));
  W := V;
  for m := 1 to #levelf do
    M1 := DegeneracyMapDomain(M, levelf[m][1]);
    if Dimension(M1) ne 0 then
      degeneracy_mats := [* degeneracy_map1(M1, M), DegeneracyMap(M1, M, levelf[m][1]) *];
      for d_mat in degeneracy_mats do
        V1 := VectorSpace(K1, Nrows(d_mat));
        W := W meet Kernel(hom< V -> V1 | Transpose(d_mat)>);
      end for;
    end if;
  end for;

  A1 := BasisMatrix(sub< V | W >);
  B1 := Transpose(Solution(Transpose(A1), ScalarMatrix(K1, Rank(A1), 1)));
  A2, T := EchelonForm(Transpose(B1)*M`basis_matrix);
  N := New(HilModFrm);
  N`is_new := true;
  N`Weight := M`Weight;
  N`basis_matrix := A2;
  N`basis_matrix_inv := M`basis_matrix_inv*Transpose(A1)*T^-1;
  N`QCD := M`QCD;
  N`PLD := M`PLD;
  N`weight_rep := M`weight_rep;
  N`weight_base_field := K1;
  N`DirichletCharacter := M`DirichletCharacter;
  N`Dimension := Rank(A2);

  return N;
end intrinsic;


      

/* Computing Hecke operators .*/
      
intrinsic HeckeOperators(M::HilModFrm, p::RngIntElt) -> SeqEnum
  { Given a space of Hilbert modular forms M, with non-paritious weight k, and a rational prime p, this returns
    Hecke operator T_p for an inert prime, and the sequence of Hecke operators T_p, T_(P^2)*S_Q and T_(Q^2)*S_P
    where (p) = P*Q, for a split prime. }

  k := M`Weight;
  F := BaseField(M);
  O := Integers(F);
  OH := M`QCD`QuaternionOrder;
  units := M`QCD`max_order_units;
  H := Algebra(OH);
  I := Factorisation(p*O)[1][1];
  require not IsParitious(k) : "Only computes Hecke operators for spaces with non-paritious weights!";
  require IsPrime(p) : "The second argument must be a rational prime!";
  require not IsRamified(I) : "Only implemented for non ramified primes!";

  if not assigned M`Dimension then
    _ := Dimension(M);
  end if;

  MM := M`HilModSpace;
  if IsInert(I) then
    tp := get_tps(MM, p*O);
    hecke_seq := [ hecke_operator(M, tp[<1,1>]) ];
  else
    _, u1 := IsNarrowlyPrincipal(I);
    u2 := Conjugate(u1);
    TP := get_tps(MM, u1*O)[<1,1>];
    TQ := get_tps(MM, u2*O)[<1,1>];
    assert Norm(TP[1]) eq u1;
    assert Norm(TQ[1]) eq u2;
    Tp := [ s1*s2 : s1 in TP, s2 in TQ ];
    SQTP2 := [ u2*s : s in get_tps(MM, u1^2*O)[<1,1>] ];
    SPTQ2 := [ u1*s : s in get_tps(MM, u2^2*O)[<1,1>] ];
    hecke_seq := [ hecke_operator(M, tp) : tp in [ Tp, SQTP2, SPTQ2 ] ];
  end if;

  return hecke_seq;
end intrinsic;



