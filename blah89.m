restore "Precomputations/hecke_mtrx_sessions/stupid_sesh";

k := [2,2,2];
H := HeckeCharacterGroup(N, [1,2,3]);
chi := H.0;

/*
for pp in PrimesUpTo(20, F) do
  HeckeMatrix2(Gamma, N, pp, k, chi);
end for;

HeckeMatrix2(Gamma, N, "Infinity", k, chi);
*/
