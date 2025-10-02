load "config.m";
import !"Geometry/ModFrmHil/precompute.m" : get_rids, support_for_rids;

F := QuadraticField(3);
ZF := Integers(F);

B := QuaternionAlgebra(1*ZF, InfinitePlaces(F) : Optimized);
O := MaximalOrder(B);

print "right ideal classes", #RightIdealClasses(O);

// Import necessary functions for computing right ideal classes

// Create a dummy modular forms space to use get_rids
// We need this because get_rids expects a ModFrmHil object
M := New(ModFrmHil);
M`QuaternionOrder := O;
M`Level := Factorization(2*ZF)[1][1];

// Compute right ideal classes and their left orders
print "Computing right ideal classes...";
rids := get_rids(M);
print "Number of right ideal classes:", #rids;

// Extract left orders
LOs := [I`LeftOrder : I in rids];
print "Number of left orders:", #LOs;


// Compute the ideal d the same way definite.m does
A := Algebra(O);
d := M`Level/Discriminant(A);
print "Discriminant of quaternion algebra:", Discriminant(A);
print "Level of M:", M`Level;
print "Using ideal d =", d, "for splitting map (computed as Level/Discriminant)";
print "Norm of d:", Norm(d);

// Create the splitting map
print "Creating splitting map...";
import "ModFrmHil/definite.m" : _ResidueMatrixRing;
split_map := _ResidueMatrixRing(O, d);

print "\n=== ANALYSIS OF LEFT ORDERS AND THEIR UNIT GROUPS ===\n";

// Iterate through each left order
for i := 1 to #LOs do
    LO := LOs[i];
    print "Left Order", i, "of", #LOs;
    print "Left Order:", LO;
    
    // Compute unit group
    U, unit_map := UnitGroup(LO);
    units := [Algebra(LO)! unit_map(s) : s in U];
    
    print "Order of unit group:", #U;
    print "Unit group structure:", U;
    
    // Apply splitting map to each unit
    print "Images of units under splitting map:";
    for j := 1 to #units do
        u := units[j];
        split_u := split_map(u);
        print "  Unit", j, ":", u;
        print "  Image:", split_u;
        print "";
    end for;
    
    print "----------------------------------------";
    print "";
end for;

print "=== ANALYSIS COMPLETE ===";

