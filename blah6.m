load "config.m";

SetVerbose("HilbertModularForms", 3);
SetVerbose("ModFrmHil", 3);

F := QuadraticField(2);
M := GradedRingOfHMFs(F, 100);
k := [1, 2];

chi_label := "-2.0.1_196.1_2u0.1u1.2u";
chi := FullChiLabelToHeckeChar(chi_label);
N := Modulus(chi);
Mk := HMFSpace(M, N, k, chi);
Sk := CuspFormBasis(Mk : StableOnly:=true);

