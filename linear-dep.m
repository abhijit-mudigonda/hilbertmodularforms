function LinearDependenceForms(A)
  // a list of forms
  B := 100;
  P := PrimesUpTo(B);
  // matrix with rows coefficients of f
  M := Matrix(Rationals(), [[Coefficient(PowerSeries(f, B), p) : p in P] : f in A]);
  return Basis(Kernel(M));
end function;

function LinearDependence(A)
  // a list of forms
  B := 100;
  P := PrimesUpTo(B);
  // matrix with rows coefficients of f
  M := Matrix(Rationals(), [[Coefficient(f, p) : p in P] : f in A]);
  return Basis(Kernel(M));
end function;


