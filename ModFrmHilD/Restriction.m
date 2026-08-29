intrinsic RestrictionToDiagonal(f::ModFrmHilDElt,M::ModFrmHilDGRng,bb::RngOrdIdl : AsCoefficients:=false) -> Any
  {Given an HMF f of parallel weight k, returns the classical modular form of weight nk and level obtained from
  restricting the component bb of the HMF to the diagonal, as a ModFrmElt. If AsCoefficients is true, instead
  returns the SeqEnum of q-expansion coefficients of the restriction, in whatever ring they naturally live,
  without coercing into a classical ModularForms space (which only supports rational coefficients) -- this
  also works when the restriction's coefficients are not rational.}
  require #SequenceToSet(Weight(f)) eq 1: "Only defined for parallel weight.";
  F := M`BaseField;
  ZF := Integers(F);
  NN := Level(f);
  N := Integers()!(Denominator(NN)^(-1)*Generator((Denominator(NN)*NN) meet Integers()));
  D := Different(ZF);
  k :=  Weight(f)[1];
  n := #Weight(f);
  fbb := f`Components[bb];
  b := Integers()!(Denominator(bb)^(-1)*Generator((Denominator(bb)*bb) meet Integers()));

  if not AsCoefficients then
    C := BaseField(F);
    R<q> := PowerSeriesRing(C);
    restriction := R!0;
    // modForms only accepts integer coefficients
    denom := 1;
  else
    raw_coeffs := AssociativeArray();
  end if;

  prec := 0;
  for j in [0 .. Precision(fbb)] do
    tracej := PositiveElementsOfTrace(bb*D^(-1),j);
    coefficient := 0;
    for nu in tracej do
      nuRed := FunDomainRep(M, nu: CheckComponent := bb);
      has_nuRed, coeffNu := IsDefined(Coefficients(fbb), nuRed);
      if not has_nuRed then
        break j;
      end if;
      coefficient +:= coeffNu;
      if not AsCoefficients then
        denom := Lcm(denom, Denominator(coeffNu));
      end if;
    end for;
    exp := j div b;
    if AsCoefficients then
      raw_coeffs[exp] := IsDefined(raw_coeffs, exp) select raw_coeffs[exp] + coefficient else coefficient;
    else
      restriction +:= coefficient*q^exp;
    end if;
    prec +:= 1;
  end for;

  if AsCoefficients then
    if #Keys(raw_coeffs) eq 0 then
      return [];
    end if;
    max_exp := Max(Keys(raw_coeffs));
    return [IsDefined(raw_coeffs, e) select raw_coeffs[e] else 0 : e in [0 .. max_exp]];
  end if;

  modForms := ModularForms(Gamma0(N),n*k);
  return modForms!(denom*(restriction +O(q^(prec))));
end intrinsic;

intrinsic PositiveElementsOfTrace(aa::RngOrdFracIdl, t::RngIntElt) -> SeqEnum[RngOrdFracIdl]
  {
    Given aa a fractional ideal and t a nonnegative integer,
    returns the totally positive elements of aa with trace t.
  }
  basis := TraceBasis(aa);
  smallest_trace := Integers()!Trace(basis[1]);
  if (t mod smallest_trace) eq 0 then
    F := NumberField(Parent(basis[1]));
    n := Degree(F);

    B := Matrix([[Evaluate(b, v) : v in InfinitePlaces(F)] : b in basis]);
    B := t * B^-1;

    // drop the first coordinate, it'll always be t / smallest_trace
    vertices := [Rationalize(Vector([v[i] : i in [2 .. n]])) : v in Rows(B)];
    assert #vertices[1] eq n-1 and #vertices eq n;
    P := Polytope(vertices);
    pts := InteriorPoints(P);
    // put t / smallest_trace back in each vector
    x := t div smallest_trace;
    return [x * basis[1] + &+[Eltseq(pt)[i] * basis[i+1] : i in [1 .. n-1]] : pt in pts];
  else
    return [];
  end if;
end intrinsic;
