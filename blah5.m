// test the code at (N, chi) which should have dihedral forms
// As of 1/1/26, it appears to run correctly - it finds a
// 2-dimensional space for each of these CM forms, and
// the eigenforms look like CM forms. 
//
// Note that only 252 has a hope of having interesting
// oldforms since for the other examples up to 272,
// the 

load "config.m";

SetVerbose("HilbertModularForms", 3);
SetVerbose("ModFrmHil", 3);

F := QuadraticField(2);
ZF := Integers(F);
M := GradedRingOfHMFs(F, 500);
k := [1, 2];

chi_labels := [
  "-2.0.1_119.2_2u1.0u1.2u",
  "-2.0.1_196.1_2u0.1u1.2u",
  "-2.0.1_224.1_2u1.0.0.1u1.2u",
  "-2.0.1_252.1_2u1.1.0u1.2u",
  "-2.0.1_272.2_2u1.1.0u1.2u",
  "-2.0.1_272.1_2u0.0.1u1.2u "
];

chi_labels := [
  "-2.0.1_476.3_2u1.0.1u1.2u"
];

codiff := Min([<Norm(nu), nu> : nu in FunDomainReps(M)[1*ZF] | not IsZero(nu)])[2];
Norm(codiff);
for chi_label in chi_labels do
  chi := FullChiLabelToHeckeChar(chi_label);
  N := Modulus(chi);
  Mk := HMFSpace(M, N, k, chi);
  Sk := CuspFormBasis(Mk : StableOnly:=true);

  eigs := Eigenbasis(Mk, Sk : P:=10);
  f := eigs[1] / Coefficient(eigs[1], 1*ZF, codiff);
  for pp in PrimesUpTo(100, F) do
    _, pi := IsNarrowlyPrincipal(pp);
    nu := pi * codiff;
    print Norm(pp), IsZero(Coefficient(f, 1*ZF, nu));
  end for;
  print "--------";
end for;
