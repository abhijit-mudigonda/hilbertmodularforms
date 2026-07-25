load "config.m";

MAX_NORM := 400;
F := QuadraticField(5);

out := Open("ideal-labels-400.txt", "w");
for N in IdealsUpTo(MAX_NORM, F) do
  fprintf out, "%o\n", LMFDBLabel(N);
end for;
delete out;
quit;
