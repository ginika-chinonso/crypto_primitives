// circom test.circom --O0 --r1cs --sym --json










pragma circom 2.1.6;

template test() {
    signal input a;
    signal input b;
    signal input c;
    signal input x;

    signal w <== x * x;
    signal y <== a * w;
    signal z <== b * x;
    signal n <== y + z + c;
}

component main{public[a,b,c,x]} = test();

// -4 * 4 = -5 ==> -x * x = -w
// -1 * 5 = -6 ==> -a * w = -y
// -2 * 4 = -7 ==> -b * x = -z
// 8 - 7 = 9 ==> n - z = newvar
// 6 + 3 = 9 ==> y + c == newvar












// template test() {
//     signal input a;
//     signal input b;
//     signal output c;

//     c <== a * b;
// }

// component main{public[a]} = test();









































// template test() {
//     signal input a; //5
//     signal input b; //3
//     signal input c; // 2
//     signal input x; //9

//     signal w <== a + b;
//     signal j <== c + x + 6;
//     signal res <== w * j + 9;
// }

// component main{public[a]} = test();

// 5 - 2 = 1 ==> w - b = a
// 6 - 4 = 8 ==> j - x = newvar1
// -5 * 6 = 9 ==> -w * j = newvar2
// 3 + const(6) = 8 ==> c + const(6) = newvar1
// -7 + const(9) = 9 ==> -res + const(9) = newvar2

// constants
// -1 => 10
// 6 => 12
// 1 => 0
// 0 => 11
// 9 => 13

































// pragma circom 2.1.6;

// include "./node_modules/circomlib/circuits/comparators.circom";

// // Create a Quadratic Equation( ax^2 + bx + c ) verifier using the below data.
// // Use comparators.circom lib to compare results if equal

// template QuadraticEquation() {
//     signal input x;     // x value
//     signal input a;     // coeffecient of x^2
//     signal input b;     // coeffecient of x 
//     signal input c;     // constant c in equation
//     signal input res;   // Expected result of the equation
//     signal output out;  // If res is correct , then return 1 , else 0 . 

//     // your code here
//     signal x_squared <== x * x;
//     signal first_term <== a * x_squared;
//     signal second_term <== b * x;
//     signal partial_sum <== first_term + second_term;

//     component equal = IsEqual();
//     equal.in[0] <== partial_sum + c;
//     equal.in[1] <== res;

//     out <== equal.out;
// }

// component main  = QuadraticEquation();

/* INPUT = {
    "x": "2",
    "a": "1",
    "b": "5",
    "c": "4",
    "res": "18"
} */


// -2 * 2 = -7 ==> -x * x = -x_squared
// -3 * 7 = -8 ==> -a * x_squared = -first_term
// -4 * 2 = -9 ==> -b * x = -second_term
// -9 - 8 = -10 ==> -second_term - first_term = -partial_sum
// -5 + 12 = 10 ==> -c + equal.in[0] = partial_sum
// -6 = -13 ==> -res = -equal.in[1]
// -11 = -1 ==> -equal.out = -out
// 15 - 13 = -12 ==> equal.isz.in - equal.in[1] = -equal.in[0]
// -14 == -11 ==> -equal.isz.out = -equal.out
// 15 * 16 = 17 ==> equal.isz.in * equal.isz.inv = new_var
// 15 * 14 = None ==> equal.isz.in * equal.isz.out == None
// -14 + 0 = 17 ==> -equal.isz.out + 1 = newvar