// level has label 2.1
// Computed with precision = 250
// generator degree bound = 14
// relation degree bound = 28
// using Standard algorithm
// component with label 1.1
// Hilbert modular variety with label 2.2.8.1-2.1-1.1-gl-0
R<[x]> := PolynomialRing(RationalField(), [ 2, 2, 4, 10 ], <"grevlexw", \[ 2, 2, 4, 10 ]>);
A := Proj(R);
eqns := [
x[1]^10 + 1606/329*x[1]^9*x[2] - 59211/329*x[1]^8*x[2]^2 + 390216/329*x[1]^7*x[2]^3 - 1194654/329*x[1]^6*x[2]^4 + 2003940/329*x[1]^5*x[2]^5 - 1884126/329*x[1]^4*x[2]^6 + 874056/329*x[1]^3*x[2]^7 - 46251/329*x[1]^2*x[2]^8 - 120506/329*x[1]*x[2]^9 + 4943/47*x[2]^10 + 886464/329*x[1]^8*x[3] - 6787584/329*x[1]^7*x[2]*x[3] + 33972480/329*x[1]^6*x[2]^2*x[3] - 158713344/329*x[1]^5*x[2]^3*x[3] + 459492480/329*x[1]^4*x[2]^4*x[3] - 742362624/329*x[1]^3*x[2]^5*x[3] + 665397504/329*x[1]^2*x[2]^6*x[3] - 311523840/329*x[1]*x[2]^7*x[3] + 59638464/329*x[2]^8*x[3] + 2495038464/329*x[1]^6*x[3]^2 - 6063040512/47*x[1]^5*x[2]*x[3]^2 + 288487111680/329*x[1]^4*x[2]^2*x[3]^2 - 926203465728/329*x[1]^3*x[2]^3*x[3]^2 + 1434674985984/329*x[1]^2*x[2]^4*x[3]^2 - 1047451834368/329*x[1]*x[2]^5*x[3]^2 + 290439447552/329*x[2]^6*x[3]^2 + 4049854267392/329*x[1]^4*x[3]^3 - 47442460409856/329*x[1]^3*x[2]*x[3]^3 + 177909226536960/329*x[1]^2*x[2]^2*x[3]^3 - 229690488913920/329*x[1]*x[2]^3*x[3]^3 + 95173868519424/329*x[2]^4*x[3]^3 + 2148670132715520/329*x[1]^2*x[3]^4 - 10891291885830144/329*x[1]*x[2]*x[3]^4 + 9757075848560640/329*x[2]^2*x[3]^4 - 1847328768/329*x[1]^5*x[4] + 12294291456/329*x[1]^4*x[2]*x[4] - 30703878144/329*x[1]^3*x[2]^2*x[4] + 36819173376/329*x[1]^2*x[2]^3*x[4] - 21467234304/329*x[1]*x[2]^4*x[4] + 700710912/47*x[2]^5*x[4] - 2916995825664/329*x[1]^3*x[3]*x[4] + 21079422664704/329*x[1]^2*x[2]*x[3]*x[4] - 54542318174208/329*x[1]*x[2]^2*x[3]*x[4] + 36379891335168/329*x[2]^3*x[3]*x[4] - 169075682574336/47*x[1]*x[3]^2*x[4] + 3212437968912384/329*x[2]*x[3]^2*x[4] + 1014454095446016/329*x[4]^2
];
S := Scheme(A,eqns);
// Computation took 725.560 seconds
// Sanity check passed: Hilbert series agree!
// Total computation for all components took 730.190 seconds

