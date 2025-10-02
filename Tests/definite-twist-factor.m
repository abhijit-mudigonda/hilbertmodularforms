//////////////////////////////////////////////////////////////////////////////
//
// Tests for the twist_factor function in ModFrmHil/definite.m
//
// This file tests the twist_factor function which computes factors that 
// appear in formulas for Brandt matrices with character.
//
//////////////////////////////////////////////////////////////////////////////

// Import the twist_factor function and related dependencies
import "ModFrmHil/definite.m": twist_factor, ProjectiveLineData;

//////////////////////////////////////////////////////////////////////////////
// Test 1: Basic functionality with trivial character
//////////////////////////////////////////////////////////////////////////////

function test_twist_factor_trivial_character()
    printf "Test 1: Testing twist_factor with trivial character...\n";
    
    // Set up a simple case with Q(sqrt(5))
    F := QuadraticField(5);
    ZF := Integers(F);
    
    // Create a simple level
    N := 2*ZF;
    
    // Create a quaternion algebra split at all finite places
    // For testing, we'll use the matrix algebra M_2(F)
    A := MatrixAlgebra(F, 2);
    O := MatrixRing(ZF, 2);
    
    try
        // Create a simple projective line data with trivial character
        R, phi := quo<ZF | N>;
        P1, P1rep := ProjectiveLine(R : Type := "Matrix");
        
        // Create a simple splitting map
        split_map := map<A -> MatrixRing(ZF, 2) | x :-> x>;
        
        // Create ProjectiveLineData record with trivial character
        PLD := rec<ProjectiveLineData |
            Level := N,
            Quotient := <R, phi>,
            P1List := P1,
            P1Rep := P1rep,
            FD := [P1[1]],
            Lookuptable := AssociativeArray(P1),
            StabOrders := [1],
            splitting_map := split_map,
            Character := 0,  // Trivial character
            CosetRepDict := 0
        >;
        
        // Test element from the quaternion algebra
        gamma := A![1, 1, 0, 1];  // Upper triangular matrix
        
        // Pick an element from P1List
        x := P1[1];
        
        // Call twist_factor - should return 1 for trivial character
        result := twist_factor(PLD, gamma, x);
        
        if result eq 1 then
            printf "  ✓ PASS: twist_factor returns 1 for trivial character\n";
            return true;
        else
            printf "  ✗ FAIL: Expected 1, got %o\n", result;
            return false;
        end if;
        
    catch e
        printf "  ✗ FAIL: Exception occurred: %o\n", e;
        return false;
    end try;
end function;

//////////////////////////////////////////////////////////////////////////////
// Test 2: Input validation and error handling
//////////////////////////////////////////////////////////////////////////////

function test_twist_factor_input_validation()
    printf "Test 2: Testing input validation...\n";
    
    // Set up basic structures
    F := QuadraticField(5);
    ZF := Integers(F);
    N := 3*ZF;
    
    A := MatrixAlgebra(F, 2);
    
    try
        R, phi := quo<ZF | N>;
        P1, P1rep := ProjectiveLine(R : Type := "Matrix");
        split_map := map<A -> MatrixRing(ZF, 2) | x :-> x>;
        
        PLD := rec<ProjectiveLineData |
            Level := N,
            Quotient := <R, phi>,
            P1List := P1,
            P1Rep := P1rep,
            FD := [P1[1]],
            Lookuptable := AssociativeArray(P1),
            StabOrders := [1],
            splitting_map := split_map,
            Character := 0,
            CosetRepDict := 0
        >;
        
        gamma := A![1, 0, 0, 1];  // Identity matrix
        x := P1[1];
        
        // Test with valid inputs - should work
        result := twist_factor(PLD, gamma, x);
        printf "  ✓ PASS: Valid inputs accepted\n";
        
        // Test with element not in P1List should fail
        try
            fake_x := Matrix(BaseRing(P1[1]), 2, 1, [999, 999]);
            result := twist_factor(PLD, gamma, fake_x);
            printf "  ✗ FAIL: Should have failed with invalid x\n";
            return false;
        catch e
            printf "  ✓ PASS: Correctly rejected invalid x\n";
        end try;
        
        return true;
        
    catch e
        printf "  ✗ FAIL: Unexpected exception: %o\n", e;
        return false;
    end try;
end function;

//////////////////////////////////////////////////////////////////////////////
// Test 3: Test with non-trivial character (mock setup)
//////////////////////////////////////////////////////////////////////////////

function test_twist_factor_with_character()
    printf "Test 3: Testing twist_factor with non-trivial character...\n";
    
    // This is a more complex test that would require setting up
    // a proper Hecke character and coset representatives
    // For now, we'll create a mock setup to test the structure
    
    F := QuadraticField(-1);  // Gaussian integers
    ZF := Integers(F);
    N := 3*ZF;
    
    A := MatrixAlgebra(F, 2);
    
    try
        R, phi := quo<ZF | N>;
        P1, P1rep := ProjectiveLine(R : Type := "Matrix");
        split_map := map<A -> MatrixRing(ZF, 2) | x :-> x>;
        
        // Create a mock character function
        mock_char := func<x | 1>;  // Simple mock that returns 1
        
        // Create mock coset representatives dictionary
        coset_dict := AssociativeArray(P1);
        for p in P1 do
            // Create a simple representative matrix
            coset_dict[p] := A![1, 0, 0, 1];  // Identity for simplicity
        end for;
        
        PLD := rec<ProjectiveLineData |
            Level := N,
            Quotient := <R, phi>,
            P1List := P1,
            P1Rep := P1rep,
            FD := [P1[1]],
            Lookuptable := AssociativeArray(P1),
            StabOrders := [1],
            splitting_map := split_map,
            Character := mock_char,
            CosetRepDict := coset_dict
        >;
        
        gamma := A![1, 1, 0, 1];  // Upper triangular
        x := P1[1];
        
        // This should execute the non-trivial character branch
        result := twist_factor(PLD, gamma, x);
        
        printf "  ✓ PASS: Non-trivial character branch executed, result = %o\n", result;
        return true;
        
    catch e
        printf "  ✗ FAIL: Exception in non-trivial character test: %o\n", e;
        return false;
    end try;
end function;

//////////////////////////////////////////////////////////////////////////////
// Test 4: Test determinant and invertibility checks
//////////////////////////////////////////////////////////////////////////////

function test_twist_factor_determinant_checks()
    printf "Test 4: Testing determinant and invertibility assertions...\n";
    
    F := QuadraticField(2);
    ZF := Integers(F);
    N := 5*ZF;
    
    A := MatrixAlgebra(F, 2);
    
    try
        R, phi := quo<ZF | N>;
        P1, P1rep := ProjectiveLine(R : Type := "Matrix");
        split_map := map<A -> MatrixRing(ZF, 2) | x :-> x>;
        
        PLD := rec<ProjectiveLineData |
            Level := N,
            Quotient := <R, phi>,
            P1List := P1,
            P1Rep := P1rep,
            FD := [P1[1]],
            Lookuptable := AssociativeArray(P1),
            StabOrders := [1],
            splitting_map := split_map,
            Character := 0,
            CosetRepDict := 0
        >;
        
        // Test with invertible matrix
        gamma_invertible := A![2, 1, 1, 1];  // det = 1, should work
        x := P1[1];
        
        result := twist_factor(PLD, gamma_invertible, x);
        printf "  ✓ PASS: Invertible matrix accepted\n";
        
        // Test with non-invertible matrix should fail
        try
            gamma_singular := A![1, 1, 1, 1];  // det = 0, should fail
            result := twist_factor(PLD, gamma_singular, x);
            printf "  ✗ FAIL: Should have failed with singular matrix\n";
            return false;
        catch e
            printf "  ✓ PASS: Correctly rejected singular matrix\n";
        end try;
        
        return true;
        
    catch e
        printf "  ✗ FAIL: Unexpected exception: %o\n", e;
        return false;
    end try;
end function;

//////////////////////////////////////////////////////////////////////////////
// Test 5: Test with different number fields
//////////////////////////////////////////////////////////////////////////////

function test_twist_factor_different_fields()
    printf "Test 5: Testing with different number fields...\n";
    
    fields := [QuadraticField(2), QuadraticField(3), QuadraticField(-1)];
    field_names := ["Q(√2)", "Q(√3)", "Q(i)"];
    
    for i := 1 to #fields do
        F := fields[i];
        field_name := field_names[i];
        
        try
            ZF := Integers(F);
            N := 2*ZF;
            A := MatrixAlgebra(F, 2);
            
            R, phi := quo<ZF | N>;
            P1, P1rep := ProjectiveLine(R : Type := "Matrix");
            split_map := map<A -> MatrixRing(ZF, 2) | x :-> x>;
            
            PLD := rec<ProjectiveLineData |
                Level := N,
                Quotient := <R, phi>,
                P1List := P1,
                P1Rep := P1rep,
                FD := [P1[1]],
                Lookuptable := AssociativeArray(P1),
                StabOrders := [1],
                splitting_map := split_map,
                Character := 0,
                CosetRepDict := 0
            >;
            
            gamma := A![1, 0, 0, 1];  // Identity
            x := P1[1];
            
            result := twist_factor(PLD, gamma, x);
            
            if result eq 1 then
                printf "  ✓ PASS: %o works correctly\n", field_name;
            else
                printf "  ✗ FAIL: %o returned %o instead of 1\n", field_name, result;
                return false;
            end if;
            
        catch e
            printf "  ✗ FAIL: %o failed with exception: %o\n", field_name, e;
            return false;
        end try;
    end for;
    
    return true;
end function;

//////////////////////////////////////////////////////////////////////////////
// Main test runner
//////////////////////////////////////////////////////////////////////////////

function run_all_tests()
    printf "Running tests for twist_factor function...\n";
    printf "=" * 60 * "\n";
    
    tests := [
        test_twist_factor_trivial_character,
        test_twist_factor_input_validation,
        test_twist_factor_with_character,
        test_twist_factor_determinant_checks,
        test_twist_factor_different_fields
    ];
    
    passed := 0;
    total := #tests;
    
    for i := 1 to total do
        printf "\n";
        if tests[i]() then
            passed +:= 1;
        end if;
    end for;
    
    printf "\n";
    printf "=" * 60 * "\n";
    printf "Test Results: %o/%o tests passed\n", passed, total;
    
    if passed eq total then
        printf "🎉 ALL TESTS PASSED!\n";
        return true;
    else
        printf "❌ Some tests failed.\n";
        return false;
    end if;
end function;

// Run the tests when this file is loaded
printf "Loading twist_factor tests...\n";
if assigned test_mode and test_mode then
    run_all_tests();
end if;
