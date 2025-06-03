load "config.m";
AttachSpec("~/Parallel.magma/spec");
import "ModFrmHil/diamond.m" : SaveHeckeMatrix; 
// SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

B := QuaternionAlgebra<F | -3, -12*F.1>;
O := MaximalOrder(B);
Gamma := FuchsianGroup(O);

/*
N := 3*ZF;
_ := cIdealDatum(Gamma, N);

k := [2,2,4];
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.0;
*/

N := ideal<ZF | [43, 0, 0], [15, 1, 0]>;
k := [2,2,3];
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.1;
assert Order(chi) eq 2;
assert HeckeCharLabel(chi) eq "1.-2.-1.1_43.2_2u1u1.2.3u";


MIN_NORM := 1;
MAX_NORM := 111;
pps := [pp : pp in PrimesUpTo(MAX_NORM, F) | Norm(pp) ge MIN_NORM];
jobs := 8;

SaveHeckeMatrix(Gamma, N, "Infinity", k, chi);

f := function(pp)
  print "pp", IdealOneLine(pp);
  SaveHeckeMatrix(Gamma, N, pp, k, chi);
  Remove(~(Gamma`ideal_data), N);
  return 1;
end function;

function foo()
  p := ParallelCall(jobs, f, [<pp> : pp in pps], 1);
  return p;
end function;


// foo();
/*
foo();
print "==================================================================================";
print "LETS GOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO";
M := GradedRingOfHMFs(F, MAX_NORM);
Mk := HMFSpace(M, N, k, chi);
CuspFormBasis(Mk);
*/
