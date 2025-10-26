// test the code at (N, chi) which should have dihedral forms

load "config.m";

SetVerbose("HilbertModularForms", 3);
SetVerbose("ModFrmHil", 3);

F := QuadraticField(2);
M := GradedRingOfHMFs(F, 500);
k := [1, 2];

chi_labels := [
//   "-2.0.1_119.2_2u1.0u1.2u",
//  "-2.0.1_196.1_2u0.1u1.2u",
  "-2.0.1_224.1_2u1.0.0.1u1.2u",
  "-2.0.1_252.1_2u1.1.0u1.2u",
  "-2.0.1_272.2_2u1.1.0u1.2u",
  "-2.0.1_272.1_2u0.0.1u1.2u "];

for chi_label in chi_labels do
  chi := FullChiLabelToHeckeChar(chi_label);
  N := Modulus(chi);
  Mk := HMFSpace(M, N, k, chi);
  Sk := CuspFormBasis(Mk : StableOnly:=true);
  print #Sk;
end for;
