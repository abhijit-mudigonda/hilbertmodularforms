load "config.m";

SetVerbose("ModFrmHil", 3);
SetVerbose("HilbertModularForms", 3);

// chi_label := "-2.0.1_32.1_2u0.1.1u1.2u";
chi_label := "-2.0.1_98.2_2u1.1u1.2u";
chi_label := "-29.0.1_91.4_2u1.0u1.2u";


k := [1, 2];
chi := FullChiLabelToHeckeChar(chi_label);
print "chi", chi;
N := Modulus(chi);
print "N", IdealOneLine(N);
print "conductor", IdealOneLine(Conductor(chi));
F := NumberField(Order(N));

assert IsCompatibleWeight(chi, k);

M := GradedRingOfHMFs(F, 350);
Mk := HMFSpace(M, N, k, chi);
Sk := CuspFormBasis(Mk : StableOnly:=true);

