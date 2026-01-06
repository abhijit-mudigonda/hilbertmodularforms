// CAN DELETE
load "config.m";

SetVerbose("HilbertModularForms", 3);
SetVerbose("ModFrmHil", 3);

F := QuadraticField(2);
M := GradedRingOfHMFs(F, 500);
k := [2, 3];

chi_labels := [
  "-2.0.1_119.2_2u1.0u1.2u",
  "-2.0.1_196.1_2u0.1u1.2u",
  "-2.0.1_224.1_2u1.0.0.1u1.2u",
  "-2.0.1_252.1_2u1.1.0u1.2u",
  "-2.0.1_272.2_2u1.1.0u1.2u",
  "-2.0.1_272.1_2u0.0.1u1.2u",
  "-2.0.1_391.1_2u1.0u1.2u",
  "-2.0.1_441.1_2u0.1u1.2u",
  "-2.0.1_476.3_2u1.0.1u1.2u",
  "-2.0.1_511.2_2u0.1u1.2u",
  "-2.0.1_512.1_2u0.1.1u1.2u",
  "-2.0.1_527.1_2u0.1u1.2u",
  "-2.0.1_544.2_2u0.0.1.1u1.2u",
  "-2.0.1_544.1_2u1.1.1.0u1.2u"
];

for chi_label in chi_labels do
  chi_12 := FullChiLabelToHeckeChar(chi_label);
  H := Parent(chi_12);
  N := Modulus(chi_12);
  chis := [];
  for chi in Elements(H) do
    if IsCompatibleWeight(chi, k) then
      Append(~chis, chi);
    end if;
  end for;

  for chi in chis do
    for tup in Factorization(N) do
      if tup[2] eq 1 and ((N / tup[1]) subset Conductor(chi)) then
        print IdealOneLine(tup[1]), HeckeCharLabel(chi);
      end if;
    end for;
  end for;
end for;
