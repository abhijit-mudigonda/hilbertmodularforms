
///// Shintani Algorithms + Enumerations of Totally positive elements in ideals /////////
// Todo Massive clean up

// Helper Functions 
intrinsic EmbedNumberField(nu::RngOrdElt, places::SeqEnum) -> SeqEnum
  { Input: nu an element of ZF where F is totally real
    Output: A tuple of the real embeddings of nu in RR}
  return [Evaluate(nu, pl) : pl in places];
end intrinsic;


intrinsic Slope(alpha::RngOrdElt) -> FldReElt
  { Input:  alpha, an element of ZF for F a totally real quadratic number field
    Output: The "slope" defined by alpha: sigma_2(alpha)/sigma_1(alpha) where sigma_i is the ith embedding of F}
  OK := Parent(alpha);
  K := NumberField(OK);
  places := InfinitePlaces(K);
  return Evaluate(alpha, places[2]) / Evaluate(alpha, places[1]);
end intrinsic;


// Rearranges the basis for an ideal so that the second basis vector has trace 0
intrinsic TraceBasis(bb::RngOrdFracIdl) -> SeqEnum
  {Given a fractional ideal bb, returns a basis (a,b) in Smith normal form where Trace(a) = n and Trace(b) = 0}
  basis := Basis(bb);
  ZF := Parent(basis[2]);
  Tr := Matrix([[Trace(basis[i]) : i in [1..#basis]]]);
  _,_,Q := SmithForm(Tr);
  ChangeofBasisMatrix := ChangeRing(Q,ZF);
  NewBasis := Eltseq(Vector(basis)*ChangeofBasisMatrix);
  return NewBasis;
end intrinsic;




///////////////////////////////////////////////////
//                                               //
//   Enumeration of Totally Positive elements    //
//                                               //
///////////////////////////////////////////////////

// Totally Positive Elements in an Ideal 
/* Idea: I've hopefully obtained basis {a,b} for the ideal bb where Tr(a) = n and Tr(b) = 0. Elements in ideal will look like xa+yb where x,y \in Z and have embedding xa_1+ yb_1 and xa_2+ yb_2.  All totally positive elements of given trace t will satisfy 
1).    t = Tr(xa+yb)    <=>   t = xn    
2).    xa+yb >> 0.     <=>   y > -x*a_1/b_1   and  y > -x*a_2/b_2  
Eq 1) determines the value for x while Eq 2) allows us to loop over values of y */

intrinsic PositiveElementsOfTrace(bb::RngOrdFracIdl, t::RngIntElt) -> SeqEnum[RngOrdFracIdl]
  {Given bb a fractional ideal, t a trace bound, returns the totally positive elements of bb with trace t.}
  Basis := TraceBasis(bb);
  places := InfinitePlaces(NumberField(Parent(Basis[1])));
  SmallestTrace := Trace(Basis[1]);
  T := [];
  if t mod SmallestTrace eq 0 then
    x := t div SmallestTrace;
    a_1 := Evaluate(Basis[1],places[1]); b_1 := Evaluate(Basis[2],places[1]);
    a_2 := Evaluate(Basis[1],places[2]); b_2 := Evaluate(Basis[2],places[2]);
    Lower := Ceiling(Min(-x*a_1/b_1,-x*a_2/b_2));
    Upper := Floor(Max(-x*a_1/b_1,-x*a_2/b_2));
    for y in [Lower .. Upper] do
      /* Append(~T, [x*Basis[1]+y*Basis[2]]); */
      Append(~T, x*Basis[1]+y*Basis[2]);
    end for;
  end if;
  return T;
end intrinsic;


///////////////////////////////////////////////////
//                                               //
//          Shintani Domain algorithms           //
//                                               //
///////////////////////////////////////////////////

// Helper Functions
// Returns the slopes of the upper and lower walls for the Shintani Domain
intrinsic ShintaniWalls(ZF::RngOrd) -> Any
  {returns lower and upper walls of the Shintani domain}
  F := NumberField(ZF);
  assert Degree(F) eq 2;
  _, F := IsQuadratic(F);
  places := InfinitePlaces(F);
  eps := FundamentalUnit(F);
  if Norm(eps) eq -1 then
    // In this case CK = CK^+ so the totally positive units are squares i.e. the subgroup generated by eps^2
    eps := eps^2;
  else 
    if not IsTotallyPositive(eps) then
      // In this case CK not equal to CK^+ so there are no units of mixed signs. If the fundamental unit is not totally positive we multiply by -1 
      eps := -1*eps;
    end if;
  end if;
  eps1 := Evaluate(eps, places[1]);
  eps2 := Evaluate(eps, places[2]);
  // always returns smallest first
  if eps1/eps2 le eps2/eps1 then
    return Sqrt(eps1/eps2), Sqrt(eps2/eps1);
  else
    return Sqrt(eps2/eps1), Sqrt(eps1/eps2);
  end if;
end intrinsic;


// Elements of the Shintani Domain with trace t
/* Idea: I've hopefully obtained basis {a,b} for the ideal bb where Tr(a) = n and Tr(b) = 0. Elements in ideal will look like xa+yb where x,y \in Z and have embedding xa_1+ yb_1 and xa_2+ yb_2. All totally positive elements of given trace t will satisfy
1).    t = Tr(xa+yb)    <=>   t = xn
2).    C_1 < (xa_1+yb_1)/(xa_2+yb_2) < C_2.     <=>   (C_1*x*a_2 -x*a_1)/(b_1-C_1*b_2) < y   and  y < (C_2*x*a_2 -x*a_1)/(b_1-C_2*b_2)
where C1 and C2 are the slope bounds on the shintani domain. Eq 1) determines the value for x while Eq 2) allows us to loop over values of y */
intrinsic ShintaniDomainOfTrace(bb::RngOrdFracIdl, t::RngIntElt) -> SeqEnum[RngOrdFracIdl]
  {Given bb a fractional ideal, t a trace bound, returns the totally positive elements of bb in the balanced Shintani cone with trace t.}  
  Basis := TraceBasis(bb);
  F := NumberField(Parent(Basis[1]));
  ZF := Integers(F);
  places := InfinitePlaces(F);
  if t eq 0 then
    return [ZF!0];
  else  
    // Orienting basis
    if Evaluate(Basis[2],places[2]) lt 0 then
      Basis := [Basis[1], -Basis[2]];
    end if;
    SmallestTrace := Trace(Basis[1]);
    T := [];
    if t mod SmallestTrace eq 0 then
      x := t div SmallestTrace;
      C1,C2 := ShintaniWalls(ZF);
      a1 := Evaluate(Basis[1],places[1]); b1 := Evaluate(Basis[2],places[1]);
      a2 := Evaluate(Basis[1],places[2]); b2 := Evaluate(Basis[2],places[2]);
      Lower := (C2*x*a2 -x*a1)/(b1-C2*b2); Upper := (C1*x*a2 -x*a1)/(b1-C1*b2); 
      // Magma has some extreme problems with .999999999 /= 1. That is why this is defined in a terrible manner. It removes points that lie on the upper wall
      prec := Precision(Lower);
      if Abs(Round(Lower) - Lower) lt 10^(-prec/2) then Lower := Round(Lower); else Lower := Ceiling(Lower); end if;
      if Abs(Round(Upper) - Upper) lt 10^(-prec/2) then Upper := Round(Upper)-1; else Upper := Floor(Upper); end if;
      for y in [Lower .. Upper] do
        Append(~T, x*Basis[1]+y*Basis[2]);
      end for;
    end if;
    return T;
  end if;
end intrinsic;
  


///////////////////////////////////////////////////
//                                               //
//          Shintani Reduction Algorithms        //
//                                               //
///////////////////////////////////////////////////


// Shintani Reduction Algorithm 1 (Currently in use)
// The Shintani Domain above is stored in an array and this looks up the ideal 
intrinsic ReduceShintani(nu::RngOrdElt, bb::RngOrdFracIdl, M::ModFrmHilD) -> SeqEnum
  {Speed up for Reduce Shintani}
  ZF := Integers(M); 
  I := nu*ZF;
  ShintaniRep := ReduceIdealToShintaniRep(M)[bb][I];
  return ShintaniRep;
end intrinsic;


// Shintani Reduction Algorithm 2 
intrinsic ReduceShintaniMinimizeTrace(nu::RngOrdElt) -> Any
  {}
  if nu eq 0 then
    return Parent(nu)!0;
  end if;
  assert IsTotallyPositive(nu);
  ZF := Parent(nu);
  F := NumberField(ZF);
  places := InfinitePlaces(F);
  eps := FundamentalUnit(ZF);
  // determine signs of eps and make eps totally positive
  eps_RR := EmbedNumberField(eps, places);
  assert #eps_RR eq 2; // only for quadratic fields right now
  pos_count := 0;
  for i := 1 to #places do
    if eps_RR[i] gt 0 then
      pos_count +:= 1;
    end if;
  end for;
  if pos_count eq 0 then
    eps := -eps;
  elif pos_count eq 1 then
    eps := eps^2;
  else
    eps := eps;
  end if;
  eps_RR := EmbedNumberField(eps, places);
  slope_eps := Slope(eps);
  slope_nu := Slope(nu);
  // TODO: do we know calculus?
  // r := -Floor( RealField(100)!(Log(slope_nu)/Log(slope_eps)) ); // old formula
  /* r := Integers()!((1/2)*Round(Log(RealField(100)!slope_nu)/Log(RealField(100)!eps_RR[1]))); */
  RR := RealField(100);
  ratio := RR!(1/2)*Log(RR!slope_nu)/Log(RR!eps_RR[1]);
  ratio_ceiling := Ceiling(ratio);
  ratio_floor := Floor(ratio);
  nu_ceiling := eps^ratio_ceiling*nu;
  nu_floor := eps^ratio_floor*nu;
  slope_nu_ceiling := Slope(nu_ceiling);
  slope_nu_floor := Slope(nu_floor);
  slopes := [slope_nu_floor, slope_nu_ceiling];
  nus := [nu_floor, nu_ceiling];
  ParallelSort(~slopes, ~nus);
  if IsShintaniReduced(nus[1]) then
    return nus[1];
  else
    assert IsShintaniReduced(nus[2]);
    return nus[2];
  end if;
end intrinsic;


// Test if an element is Shintani reduced 
intrinsic IsShintaniReduced(nu::RngOrdElt) -> BoolElt
  {}
  // zero is Shintani Reduced
  if nu eq Parent(nu)!0 then
    return true;
  end if;
  // wall1<wall2
  wall1, wall2 := ShintaniWalls(Parent(nu));
  slope := Slope(nu);
  prec := Precision(Parent(slope));
  // walls with fuzz
  if (wall1-10^(-prec/2) le slope) and (slope le wall2+10^(-prec/2)) then
    return true;
  else
    return false;
  end if;
end intrinsic;


// Conversion : Shintani elements < = > Ideals 
// Converts pairs (bb,nu) <-> (bb,n) based on the set of representatives bb for Cl^+(F) 

intrinsic IdealToShintaniRepresentative(M::ModFrmHilD, bb::RngOrdIdl, n::RngOrdIdl) -> ModFrmHilDElt
  {Takes a representative [bb] in Cl^+(F) and an integral ideal n in ZF with [n] = [bb^(-1)] and returns Shintani representative (nu) = n*bb}
  _,gen := IsPrincipal(n*bb);
  ShintaniGenerator := ReduceShintani(gen,bb);
  return ShintaniGenerator;
end intrinsic;


intrinsic ShintaniRepresentativeToIdeal(bb::RngOrdFracIdl, nu::RngOrdElt) -> ModFrmHilDElt
  {Takes a representative [bb^(-1)] in Cl^+(F) and a nu in bb_+ and returns the integral ideal n = bb^(-1)*(nu) in ZF}
  ZF := Parent(nu);
  n := nu*bb^(-1);
  return NicefyIdeal(n);
end intrinsic;




/// Assorted Shintani Shenanigans ////////

/* 
intrinsic ReduceShintaniComputeIdeal(nu::RngOrdElt, bb::RngOrdFracIdl, shintani_reps::Assoc) -> Any
  {}
  reps := [];
  for t in Keys(shintani_reps[bb]) do
    reps cat:= shintani_reps[bb][t];
  end for;
  return ReduceShintaniComputeIdeal(nu, reps);
end intrinsic;

intrinsic ReduceShintaniComputeIdeal(nu::RngOrdElt, shintani_reps::SeqEnum[RngOrdElt]) -> Any
  {}
  if nu eq 0 then
    return Parent(nu)!0;
  end if;
  assert IsTotallyPositive(nu);
  ZF := Parent(nu);
  nu_ideal := ideal<ZF|nu>;
  shintani_ideals := [ideal<ZF|a> : a in shintani_reps];
  matches := [];
  for i := 1 to #Keys(shintani_reps) do
    I := ideal<ZF|shintani_reps[i]>;
    if nu_ideal eq I then
      Append(~matches, [* I, i *]);
    end if;
  end for;
  assert #matches eq 1;
  return shintani_reps[matches[1][2]];
end intrinsic;
*/

/* intrinsic ReduceShintani(nu::RngOrdElt, bb::RngOrdFracIdl, shintani_reps::Assoc) -> Any */
  /* {} */
  /* nu_reduced_by_ideal := ReduceShintaniComputeIdeal(nu, bb, shintani_reps); */
  /* nu_reduced_by_trace := ReduceShintaniMinimizeTrace(nu); */
  /* return nu_reduced_by_ideal; */
  // sanity check using trace when ready
  /* if nu_reduced_by_ideal eq nu_reduced_by_trace then */
  /*   return nu_reduced_by_ideal; */
  /* else */
  /*   error_message := Sprintf("Shintani walls?\n"); */
  /*   error_message *:= Sprintf("nu using ideals = %o\n", nu_reduced_by_ideal); */
  /*   error_message *:= Sprintf("Trace(nu using ideals) = %o\n", Trace(nu_reduced_by_ideal)); */
  /*   error_message *:= Sprintf("nu using trace = %o\n", nu_reduced_by_trace); */
  /*   error_message *:= Sprintf("Trace(nu using trace) = %o\n", Trace(nu_reduced_by_trace)); */
  /*   error error_message; */
  /* end if; */
/* end intrinsic; */

/*
intrinsic ReduceShintani(pair::SeqEnum, bb::RngOrdFracIdl, shintani_reps::Assoc) -> SeqEnum
  {}
  assert #pair eq 2;
  return [ReduceShintani(nu, bb, shintani_reps) : nu in pair];
end intrinsic;

intrinsic IdealToShintaniRepresentative(M::ModFrmHilD, nn::RngOrdIdl) -> Any
  {}
  Cl := NarrowClassGroup(M);
  mp := NarrowClassGroupMap(M);
  mp_inv := Inverse(mp);
  bb := mp(-mp_inv(nn));
  bl, nu := IsPrincipal(nn*bb); // will this always give a totally positive generator? probably not
  assert bl;
  assert IsTotallyPositive(nu);
  return ReduceShintani(nu, bb, ShintaniReps(M)), bb;
end intrinsic; */




