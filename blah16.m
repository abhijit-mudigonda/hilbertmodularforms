load "config.m";

chi := FullChiLabelToHeckeChar(chi_label);
N := Modulus(chi);

F := NumberField(Order(N));
n := Degree(F);

eis_wt := [1 : _ in [1 .. n]];
H := HeckeCharacterGroup(N, [1 .. n]);

compat_flag := false;
eis_flag := false;

for psi in Elements(H) do
  if IsCompatibleWeight(psi, eis_wt) then
    compat_flag := true;
  end if;

  if IsGamma1EisensteinWeight(psi, 1) then
    eis_flag := true;
  end if;

  if compat_flag ne eis_flag then
    print "huhhhhhhhhh ok you're wrong";
  end if;
end for;
