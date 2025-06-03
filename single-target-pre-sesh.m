/* 
 * This sets up a session with a common Gamma and ideal datum.
 * This session can then be loaded by every parallel branch when computing
 * Hecke matrices so that they agree on things like coset representatives,
 * representatives of the projective line, residue maps, etc.
 */

import "ModFrmHil/diamond.m" : AlgQuatEltLabel, SaveHeckeMatrix; 
load "config.m";
SetVerbose("ModFrmHil", 3);

R<x> := PolynomialRing(Rationals());
F := NumberField(x^3-x^2-2*x+1);
ZF := Integers(F);

B := QuaternionAlgebra<F | -3, -12*F.1>;
O := MaximalOrder(B);
Gamma := FuchsianGroup(O);

N := 3*ZF;
_ := cIdealDatum(Gamma, N);

/*
F_label := HackyFieldLabel(F);
alg_label := Join([AlgQuatEltLabel(B.i^2) : i in [1 .. 2]], "_");
ord_label := Join([AlgQuatEltLabel(O.i) : i in [1 .. 3]], "_");
N_label := LMFDBLabel(N);
chi_label := HeckeCharLabel(chi : full_label:=false);
*/

// apparently the save and restore commands can't take variable string names
// for some stupid reason
// 
// hms_dir := "Precomputations/hecke_mtrx_sessions/";
// filename := Join([F_label, alg_label, ord_label, N_label, chi_label], "=") cat "_hecke_mtrx_sesh";
// filepath := hms_dir cat filename;

save "Precomputations/hecke_mtrx_sessions/stupid_sesh";
