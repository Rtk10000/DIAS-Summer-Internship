#include <iostream>
#include <random>
#include <ctime>
#include <string>
#include <complex>
#include <fstream>
#include <sstream>
#include "eigen/Eigen/Dense"
#include <iomanip>

const double g = 1;
const int N = 3;
const double dt = 1e-4;
const int iterations = 1e6;

typedef std::complex<double> complex;
typedef Eigen:: Matrix<complex, N, N> matrix;
typedef Eigen:: Matrix<complex, 1, N> row;
typedef Eigen:: Matrix<complex, N, 1> col;

double colabsqr(col a)
{
    complex b = (a.adjoint()*a);
    return b.real();
}

double rowabsqr(row a)
{
    complex b = (a*a.adjoint());
    return b.real();
}

double modsqr(complex a)
{
    return pow(abs(a),2);
}

matrix commutator(matrix A, matrix B)
{
    return A*B - B*A;
}

matrix anticommutator(matrix A, matrix B)
{
    return A*B + B*A;
}

matrix gamma(col Z41, col Z42, col Z43, row Z14, row Z24, row Z34, matrix phi41, matrix phi42, matrix phi43, double c4)
{
    return Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34
    + commutator(phi41,phi41.adjoint()) + commutator(phi42,phi42.adjoint()) + commutator(phi43,phi43.adjoint()) - c4/(g*g) * matrix::Identity();
} 

double K(
    complex U12, complex U13, complex U21, complex U23, complex U31, complex U32, row U14, row U24, row U34, col U41, col U42, col U43,
    double V11, double V12, double V13, double V14, double V15, double V16, double V17, double V18, double V19,
    double V21, double V22, double V23, double V24, double V25, double V26, double V27, double V28, double V29,
    double V31, double V32, double V33, double V34, double V35, double V36, double V37, double V38, double V39,
    matrix V41, matrix V42, matrix V43, matrix V44, matrix V45, matrix V46, matrix V47, matrix V48, matrix V49
    )
{
    // The piece pertaining to the the Z velocities (U's)
    double scalar_sum = modsqr(U12) + modsqr(U13) + modsqr(U21) + modsqr(U23) + modsqr(U31) + modsqr(U32);
    double row_sum = (U14.adjoint()*U14 + U24.adjoint()*U24 + U34.adjoint()*U34).trace().real();
    double col_sum = (U41.adjoint()*U41 + U42.adjoint()*U42 + U43.adjoint()*U43).trace().real();

    // The piece pertaining to the X velocities (V's)
    
    // for k=1
    double sum1 = pow(V11,2) + pow(V12,2) + pow(V13,2) + pow(V14,2) + pow(V15,2) + pow(V16,2) + pow(V17,2) + pow(V18,2) + pow(V19,2);
    // for k=2 
    double sum2 = pow(V21,2) + pow(V22,2) + pow(V23,2) + pow(V24,2) + pow(V25,2) + pow(V26,2) + pow(V27,2) + pow(V28,2) + pow(V29,2);
    // for k=3
    double sum3 = pow(V31,2) + pow(V32,2) + pow(V33,2) + pow(V34,2) + pow(V35,2) + pow(V36,2) + pow(V37,2) + pow(V38,2) + pow(V39,2);
    // for k=4
    double sum4 = (V41*V41 + V42*V42 + V43*V43 + V44*V44 + V45*V45 + V46*V46 + V47*V47 + V48*V48 + V49*V49).trace().real();

    return 0.5 * (scalar_sum + row_sum + col_sum + sum1 + sum2 + sum3 + sum4);
}

//Need to add V_gauge to the potential

double V_D(
    complex Z12, complex Z13, complex Z21, complex Z23, complex Z31, complex Z32, row Z14, row Z24, row Z34, col Z41, col Z42, col Z43,
    double c1, double c2, double c3, double c4,
    matrix phi41, matrix phi42, matrix phi43)
{
    // Starting with V_D {explicitly summing each k gives 4 sums}

    // For k = 1
    double sum1 = pow(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) - c1/g*g, 2);
    // For k = 2
    double sum2 = pow(modsqr(Z21) + modsqr(Z23) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/g*g, 2);
    // For k = 3
    double sum3 = pow(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(Z13) - modsqr(Z23) - colabsqr(Z43) - c3/g*g, 2);
    // For k = 4
    double sum4 = ((
        Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34 +
        commutator(phi41,phi41.adjoint()) + commutator(phi42,phi42.adjoint()) + commutator(phi43,phi43.adjoint()) - (c4/g*g) * matrix::Identity()
        )*(
        Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34 +
        commutator(phi41,phi41.adjoint()) + commutator(phi42,phi42.adjoint()) + commutator(phi43,phi43.adjoint()) - (c4/g*g) * matrix::Identity()
        )).trace().real();    
    
    return 0.5 * (sum1 + sum2 + sum3 + sum4);
}

double V_gauge(
    complex Z12, complex Z13, complex Z21, complex Z23, complex Z31, complex Z32, row Z14, row Z24, row Z34, col Z41, col Z42, col Z43,
    double X11, double X12, double X13,
    double X21, double X22, double X23,
    double X31, double X32, double X33,
    matrix X41, matrix X42, matrix X43,
    complex phi11, complex phi12, complex phi13,
    complex phi21, complex phi22, complex phi23,
    complex phi31, complex phi32, complex phi33,
    matrix phi41, matrix phi42, matrix phi43
    )
{
    // Start with the all of the definitions {Huge scope for error in indices when defining all of these variables}
    // Define H's
    complex H121 = X11 * Z12 - Z12 * X21; complex H122 = X12 * Z12 - Z12 * X22; complex H123 = X13 * Z12 - Z12 * X23;
    complex H131 = X11 * Z13 - Z13 * X31; complex H132 = X12 * Z13 - Z13 * X32; complex H133 = X13 * Z13 - Z13 * X33;
        row H141 = X11 * Z14 - Z14 * X41;     row H142 = X12 * Z14 - Z14 * X42;     row H143 = X13 * Z14 - Z14 * X43;

    complex H211 = X21 * Z21 - Z21 * X11; complex H212 = X22 * Z21 - Z21 * X12; complex H213 = X23 * Z21 - Z21 * X13;
    complex H231 = X21 * Z23 - Z23 * X31; complex H232 = X22 * Z23 - Z23 * X32; complex H233 = X23 * Z23 - Z23 * X33;
        row H241 = X21 * Z24 - Z24 * X41;     row H242 = X22 * Z24 - Z24 * X42;     row H243 = X23 * Z24 - Z24 * X43;

    complex H311 = X31 * Z31 - Z31 * X11; complex H312 = X32 * Z31 - Z31 * X12; complex H313 = X33 * Z31 - Z31 * X13;
    complex H321 = X31 * Z32 - Z32 * X21; complex H322 = X32 * Z32 - Z32 * X22; complex H323 = X33 * Z32 - Z32 * X23;
        row H341 = X31 * Z34 - Z34 * X41;     row H342 = X32 * Z34 - Z34 * X42;     row H343 = X33 * Z34 - Z34 * X43;

    col H411 = X41 * Z41 - Z41 * X11; col H412 = X42 * Z41 - Z41 * X12; col H413 = X43 * Z41 - Z41 * X13;
    col H421 = X41 * Z42 - Z42 * X21; col H422 = X42 * Z42 - Z42 * X22; col H423 = X43 * Z42 - Z42 * X23;
    col H431 = X41 * Z43 - Z43 * X31; col H432 = X42 * Z43 - Z43 * X32; col H433 = X43 * Z43 - Z43 * X33;

    // Define Y's {only non-zero for k=4}
    matrix Y411 = commutator(X41,phi41); matrix Y412 = commutator(X41,phi42); matrix Y413 = commutator(X41,phi43);
    matrix Y421 = commutator(X42,phi41); matrix Y422 = commutator(X42,phi42); matrix Y423 = commutator(X42,phi43);
    matrix Y431 = commutator(X43,phi41); matrix Y432 = commutator(X43,phi42); matrix Y433 = commutator(X43,phi43);

    //Define the X's {commutator X's with three indices}
    matrix X411 = commutator(X41,X41); matrix X412 = commutator(X41,X42); matrix X413 = commutator(X41,X43);
    matrix X421 = commutator(X42,X41); matrix X422 = commutator(X42,X42); matrix X423 = commutator(X42,X43);
    matrix X431 = commutator(X43,X41); matrix X432 = commutator(X43,X42); matrix X433 = commutator(X43,X43);
    
    double V1_gauge = (conj(H121)*H121).real() + (conj(H122)*H122).real() + (conj(H123)*H123).real() +
                      (conj(H131)*H131).real() + (conj(H132)*H132).real() + (conj(H133)*H133).real() +

                      (conj(H211)*H211).real() + (conj(H212)*H212).real() + (conj(H213)*H213).real() +
                      (conj(H231)*H231).real() + (conj(H232)*H232).real() + (conj(H233)*H233).real() +

                      (conj(H311)*H311).real() + (conj(H312)*H312).real() + (conj(H313)*H313).real() + 
                      (conj(H321)*H321).real() + (conj(H322)*H322).real() + (conj(H323)*H323).real() +

                      (H141.adjoint()*H141).trace().real() + (H142.adjoint()*H142).trace().real() + (H143.adjoint()*H143).trace().real() +
                      (H241.adjoint()*H241).trace().real() + (H242.adjoint()*H242).trace().real() + (H243.adjoint()*H243).trace().real() +
                      (H341.adjoint()*H341).trace().real() + (H342.adjoint()*H342).trace().real() + (H343.adjoint()*H343).trace().real() +
                      
                      (H411.adjoint()*H411).trace().real() + (H412.adjoint()*H412).trace().real() + (H413.adjoint()*H413).trace().real() + 
                      (H421.adjoint()*H421).trace().real() + (H422.adjoint()*H422).trace().real() + (H423.adjoint()*H423).trace().real() + 
                      (H431.adjoint()*H431).trace().real() + (H432.adjoint()*H432).trace().real() + (H433.adjoint()*H433).trace().real();

    double V2_gauge = (Y411.adjoint()*Y411).trace().real() + (Y412.adjoint()*Y412).trace().real() + (Y413.adjoint()*Y413).trace().real() + 
                      (Y421.adjoint()*Y421).trace().real() + (Y422.adjoint()*Y422).trace().real() + (Y423.adjoint()*Y423).trace().real() + 
                      (Y431.adjoint()*Y431).trace().real() + (Y432.adjoint()*Y432).trace().real() + (Y433.adjoint()*Y433).trace().real();

    double V3_gauge = (X411.adjoint()*X411).trace().real() + (X412.adjoint()*X412).trace().real() + (X413.adjoint()*X413).trace().real() +
                      (X421.adjoint()*X421).trace().real() + (X422.adjoint()*X422).trace().real() + (X423.adjoint()*X423).trace().real() +
                      (X431.adjoint()*X431).trace().real() + (X432.adjoint()*X432).trace().real() + (X433.adjoint()*X433).trace().real();
    
    return V1_gauge + V2_gauge + 0.25 * V3_gauge;
}

std:: mt19937 rng(std::time(nullptr));
std:: normal_distribution<long double> gauss_dist(0, 1);

// Initialise the Z variables only and let U be the \dot{Z} coordinates
complex Z12 = complex(gauss_dist(rng),gauss_dist(rng));
complex Z13 = complex(gauss_dist(rng),gauss_dist(rng));
complex Z21 = complex(gauss_dist(rng),gauss_dist(rng));
complex Z23 = complex(gauss_dist(rng),gauss_dist(rng)); 
complex Z31 = complex(gauss_dist(rng),gauss_dist(rng)); 
complex Z32 = complex(gauss_dist(rng),gauss_dist(rng)); 

col Z41 = col::Random();//col::Random();
col Z42 = col::Random();//col::Random();
col Z43 = col::Random();//col::Random();
row Z14 = row::Random();
row Z24 = row::Random(); 
row Z34 = row::Random();

// initialise the Z_n'

complex Z12_n = complex(0,0);
complex Z13_n = complex(0,0);
complex Z21_n = complex(0,0);
complex Z23_n = complex(0,0); 
complex Z31_n = complex(0,0); 
complex Z32_n = complex(0,0); 

col Z41_n = col::Zero();
col Z42_n = col::Zero();
col Z43_n = col::Zero();
row Z14_n = row::Zero();
row Z24_n = row::Zero();
row Z34_n = row::Zero();

// initiailise the X's
double X11 = gauss_dist(rng);
double X12 = gauss_dist(rng);
double X13 = gauss_dist(rng);
double X14 = 0;//gauss_dist(rng);
double X15 = 0;//gauss_dist(rng);
double X16 = 0;//gauss_dist(rng);
double X17 = 0;//gauss_dist(rng);
double X18 = 0;//gauss_dist(rng);
double X19 = 0;//gauss_dist(rng);

double X21 = gauss_dist(rng);
double X22 = gauss_dist(rng);
double X23 = gauss_dist(rng);
double X24 = 0;//gauss_dist(rng);
double X25 = 0;//gauss_dist(rng);
double X26 = 0;//gauss_dist(rng);
double X27 = 0;//gauss_dist(rng);
double X28 = 0;//gauss_dist(rng);
double X29 = 0;//gauss_dist(rng);

double X31 = gauss_dist(rng);
double X32 = gauss_dist(rng);
double X33 = gauss_dist(rng);
double X34 = 0;//gauss_dist(rng);
double X35 = 0;//gauss_dist(rng);
double X36 = 0;//gauss_dist(rng);
double X37 = 0;//gauss_dist(rng);
double X38 = 0;//gauss_dist(rng);
double X39 = 0;//gauss_dist(rng);

// Declaring the matrices 
matrix X41;
matrix X42; 
matrix X43; 
matrix X44; 
matrix X45; 
matrix X46; 
matrix X47;
matrix X48;
matrix X49;

// Declare the M matrices 
matrix M41 = Z41 * Z41.adjoint(); 
matrix M42 = Z42 * Z42.adjoint(); 
matrix M43 = Z43 * Z43.adjoint();
matrix M14 = Z14.adjoint() * Z14; 
matrix M24 = Z24.adjoint() * Z24; 
matrix M34 = Z34.adjoint() * Z34;

matrix M41_n = matrix::Zero();
matrix M42_n = matrix::Zero();
matrix M43_n = matrix::Zero();
matrix M14_n = matrix::Zero(); 
matrix M24_n = matrix::Zero();
matrix M34_n = matrix::Zero();

complex phi11 = X14 + complex(1,0)*X15;
complex phi12 = X16 + complex(1,0)*X17;
complex phi13 = X18 + complex(1,0)*X19;

complex phi21 = X24 + complex(1,0)*X25;
complex phi22 = X26 + complex(1,0)*X27;
complex phi23 = X28 + complex(1,0)*X29;

complex phi31 = X34 + complex(1,0)*X35;
complex phi32 = X36 + complex(1,0)*X37;
complex phi33 = X38 + complex(1,0)*X39;

// Initialise the X_n's
double X11_n = 0;//gauss_dist(rng);
double X12_n = 0;//gauss_dist(rng);
double X13_n = 0;//gauss_dist(rng);
double X14_n = 0;//gauss_dist(rng);
double X15_n = 0;//gauss_dist(rng);
double X16_n = 0;//gauss_dist(rng);
double X17_n = 0;//gauss_dist(rng);
double X18_n = 0;//gauss_dist(rng);
double X19_n = 0;//gauss_dist(rng);

double X21_n = 0;//gauss_dist(rng);
double X22_n = 0;//gauss_dist(rng);
double X23_n = 0;//gauss_dist(rng);
double X24_n = 0;//gauss_dist(rng);
double X25_n = 0;//gauss_dist(rng);
double X26_n = 0;//gauss_dist(rng);
double X27_n = 0;//gauss_dist(rng);
double X28_n = 0;//gauss_dist(rng);
double X29_n = 0;//gauss_dist(rng);

double X31_n = 0;//gauss_dist(rng);
double X32_n = 0;//gauss_dist(rng);
double X33_n = 0;//gauss_dist(rng);
double X34_n = 0;//gauss_dist(rng);
double X35_n = 0;//gauss_dist(rng);
double X36_n = 0;//gauss_dist(rng);
double X37_n = 0;//gauss_dist(rng);
double X38_n = 0;//gauss_dist(rng);
double X39_n = 0;//gauss_dist(rng);

matrix X41_n = matrix::Zero();
matrix X42_n = matrix::Zero();
matrix X43_n = matrix::Zero();
matrix X44_n = matrix::Zero();
matrix X45_n = matrix::Zero();
matrix X46_n = matrix::Zero();
matrix X47_n = matrix::Zero();
matrix X48_n = matrix::Zero();
matrix X49_n = matrix::Zero();

// Initialise U's (These are velocities for the Z coordinates)
complex U12 = complex(0,0);
complex U13 = complex(0,0);
complex U21 = complex(0,0);
complex U23 = complex(0,0);
complex U31 = complex(0,0);
complex U32 = complex(0,0);

col U41 = col::Zero();
col U42 = col::Zero();
col U43 = col::Zero();
row U14 = row::Zero();
row U24 = row::Zero();
row U34 = row::Zero();

// Initialise the V's (These are the Velocities for the X coordinates)

double V11 = 0;//gauss_dist(rng);
double V12 = 0;//gauss_dist(rng);
double V13 = 0;//gauss_dist(rng);
double V14 = 0;//gauss_dist(rng);
double V15 = 0;//gauss_dist(rng);
double V16 = 0;//gauss_dist(rng);
double V17 = 0;//gauss_dist(rng);
double V18 = 0;//gauss_dist(rng);
double V19 = 0;//gauss_dist(rng);

double V21 = 0;//gauss_dist(rng);
double V22 = 0;//gauss_dist(rng);
double V23 = 0;//gauss_dist(rng);
double V24 = 0;//gauss_dist(rng);
double V25 = 0;//gauss_dist(rng);
double V26 = 0;//gauss_dist(rng);
double V27 = 0;//gauss_dist(rng);
double V28 = 0;//gauss_dist(rng);
double V29 = 0;//gauss_dist(rng);

double V31 = 0;//gauss_dist(rng);
double V32 = 0;//gauss_dist(rng);
double V33 = 0;//gauss_dist(rng);
double V34 = 0;//gauss_dist(rng);
double V35 = 0;//gauss_dist(rng);
double V36 = 0;//gauss_dist(rng);
double V37 = 0;//gauss_dist(rng);
double V38 = 0;//gauss_dist(rng);
double V39 = 0;//gauss_dist(rng);

matrix V41 = matrix::Zero();//matrix::Random();
matrix V42 = matrix::Zero();//matrix::Random();
matrix V43 = matrix::Zero();//matrix::Random();
matrix V44 = matrix::Zero();//matrix::Random();
matrix V45 = matrix::Zero();//matrix::Random();
matrix V46 = matrix::Zero();//matrix::Random();
matrix V47 = matrix::Zero();//matrix::Random();
matrix V48 = matrix::Zero();//matrix::Random();
matrix V49 = matrix::Zero();//matrix::Random();


// Initialise c's as complex but they may be real ask Swapno if things don't work
double c1 = 0;//gauss_dist(rng);
double c2 = 0;//gauss_dist(rng);
double c3 = 0;//gauss_dist(rng);
double c4 = 0;//gauss_dist(rng);

// Need to add the piece of acceleration functions pertaining to V_gauge

// Define the Force functions (defined as in dropbox -\ddot(Zij) = FZij) {Need to add the X coordinates to these force functions}
// Z functions going first
complex FZ12(
    complex Z12, complex Z13, row Z14, complex Z21, complex Z31, col Z41,
    complex Z23, row Z24, complex Z32, col Z42,
    double X11, double X12, double X13, double X21, double X22, double X23)
{
    //V_D piece
    complex V_D = complex(2,0) * (Z12*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) - c1/(g*g))
    - Z12*(modsqr(Z21) + modsqr(Z23) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/(g*g)));

    //V_gauge piece
    complex V_gauge = complex(2,0) * Z12 * (pow(X11 - X21, 2) + pow(X12 - X22, 2) + pow(X13 - X23, 2));
    
    return V_D + V_gauge;
}

complex FZ21(
    complex Z12, complex Z13, row Z14, complex Z21, complex Z31, col Z41,
    complex Z23, row Z24, complex Z32, col Z42,
    double X11, double X12, double X13, double X21, double X22, double X23)
{
    //V_D piece
    complex V_D = complex(2,0) * (- Z21*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) - c1/(g*g))
    + Z21*(modsqr(Z21) + modsqr(Z23) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/(g*g)));

    complex V_gauge = complex(2,0) * Z21 * (pow(X21 - X11, 2) + pow(X22 - X12, 2) + pow(X23 - X13, 2));

    return V_D + V_gauge;
}

complex FZ13(
    complex Z12, complex Z13, row Z14, complex Z21, complex Z31, col Z41,
    complex Z23, row Z34, complex Z32, col Z43,
    double X11, double X12, double X13, double X31, double X32, double X33)
{
    //V_D piece
    complex V_D = complex(2,0) * (Z13*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) - c1/(g*g))
    - Z13*(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(Z13) - modsqr(Z23) - colabsqr(Z43) - c3/(g*g))); 

    complex V_gauge = complex(2,0) * Z13 * (pow(X11 - X31, 2) + pow(X12 - X32, 2) + pow(X13 - X33, 2));

    return V_D + V_gauge;
}

complex FZ31(
    complex Z12, complex Z13, row Z14, complex Z21, complex Z31, col Z41,
    complex Z23, row Z34, complex Z32, col Z43,
    double X11, double X12, double X13, double X31, double X32, double X33)
{
    //V_D piece
    complex V_D = complex(2,0) * (- Z31*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) - c1/(g*g))
    + Z31*(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(Z13) - modsqr(Z23) - colabsqr(Z43) - c3/(g*g)));

    complex V_gauge = complex(2,0) * Z31 * (pow(X31 - X11,2) + pow(X32 - X12, 2) + pow(X33 - X13, 2));
    
    return V_D + V_gauge;
}

complex FZ23(
    complex Z23, complex Z21, row Z24, complex Z12, complex Z32, col Z42,
    complex Z31, row Z34, complex Z13, col Z43,
    double X21, double X22, double X23, double X31, double X32, double X33)
{
    //V_D piece
    complex V_D = complex(2,0) * (Z23*(modsqr(Z23) + modsqr(Z21) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/(g*g))
    - Z23*(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(Z13) - modsqr(Z23) - colabsqr(Z43) - c3/(g*g)));

    complex V_gauge = complex(2,0) * Z23 * (pow(X21 - X31, 2) + pow(X22 - X32, 2) + pow(X23 - X33, 2));
    
    return V_D + V_gauge;
}

complex FZ32(
    complex Z23, complex Z21, row Z24, complex Z12, complex Z32, col Z42,
    complex Z31, row Z34, complex Z13, col Z43,
    double X21, double X22, double X23, double X31, double X32, double X33)
{
    //V_D piece
    complex V_D = complex(2,0) * (- Z32*(modsqr(Z23) + modsqr(Z21) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/(g*g))
    + Z32*(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(Z13) - modsqr(Z23) - colabsqr(Z43) - c3/(g*g)));

    complex V_gauge = complex(2,0) * Z32 * (pow(X31 - X21, 2) + pow(X32 - X22, 2) + pow(X33 - X23, 2));

    return V_D + V_gauge;
}

row FZ14(
    row Z14, complex Z12, complex Z13, complex Z21, complex Z31, col Z41, 
    col Z42, col Z43, row Z24, row Z34, matrix phi41, matrix phi42, matrix phi43,
    double X11, double X12, double X13, matrix X41, matrix X42, matrix X43)
{
    //V_D piece
    row V_D = complex(2,0) * (Z14*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) -c1/(g*g))
    - Z14*(Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34
    + commutator(phi41,phi41.adjoint()) + commutator(phi42,phi42.adjoint()) + commutator(phi43,phi43.adjoint()) - c4/(g*g) * matrix::Identity())); 

    row V_gauge = complex(2,0) * Z14 * (
        X11*X11*matrix::Identity() - 2*X11*X41 + X41*X41 +
        X12*X12*matrix::Identity() - 2*X12*X42 + X42*X42 + 
        X13*X13*matrix::Identity() - 2*X13*X43 + X43*X43
    );

    return V_D + V_gauge;
}

col FZ41(
    row Z14, complex Z12, complex Z13, complex Z21, complex Z31, col Z41, 
    col Z42, col Z43, row Z24, row Z34, matrix phi41, matrix phi42, matrix phi43,
    double X11, double X12, double X13, matrix X41, matrix X42, matrix X43)
{
    //V_D piece
    col V_D = complex(2,0) * (- Z41*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) -c1/(g*g))
    + (Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34
    + commutator(phi41,phi41.adjoint()) + commutator(phi42,phi42.adjoint()) + commutator(phi43,phi43.adjoint()) - c4/(g*g) * matrix::Identity())*Z41);

    col V_gauge = complex(2,0) * (
        X11*X11*matrix::Identity() - 2*X11*X41 + X41*X41 +
        X12*X12*matrix::Identity() - 2*X12*X42 + X42*X42 + 
        X13*X13*matrix::Identity() - 2*X13*X43 + X43*X43
    ) * Z41;

    return V_D + V_gauge;
}

row FZ24(
    row Z24, complex Z21, complex Z23, complex Z12, complex Z32, col Z41,
    col Z42, col Z43, row Z14, row Z34, matrix phi41, matrix phi42, matrix phi43,
    double X21, double X22, double X23, matrix X41, matrix X42, matrix X43)
{
    //V_D piece
    row V_D = complex(2,0) * (Z24*(modsqr(Z21) + modsqr(Z23) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/(g*g))
    - Z24*(Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34
    + commutator(phi41,phi41.adjoint()) + commutator(phi42,phi42.adjoint()) + commutator(phi43,phi43.adjoint()) - c4/(g*g) * matrix::Identity()));

    row V_gauge = complex(2,0) * Z24 * (
        X21*X21*matrix::Identity() - 2*X21*X41 + X41*X41 +
        X22*X22*matrix::Identity() - 2*X22*X42 + X42*X42 + 
        X23*X23*matrix::Identity() - 2*X23*X43 + X43*X43
    );

    return V_D + V_gauge;
}

col FZ42(
    row Z24, complex Z21, complex Z23, complex Z12, complex Z32, col Z41,
    col Z42, col Z43, row Z14, row Z34, matrix phi41, matrix phi42, matrix phi43,
    double X21, double X22, double X23, matrix X41, matrix X42, matrix X43)
{
    //V_D piece
    col V_D = complex(2,0) * (- Z42*(modsqr(Z21) + modsqr(Z23) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/(g*g))
    + (Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34
    + commutator(phi41,phi41.adjoint()) + commutator(phi42,phi42.adjoint()) + commutator(phi43,phi43.adjoint()) - c4/(g*g) * matrix::Identity())*Z42);

    col V_gauge = complex(2,0) * (
        X21*X21*matrix::Identity() - 2*X21*X41 + X41*X41 +
        X22*X22*matrix::Identity() - 2*X22*X42 + X42*X42 + 
        X23*X23*matrix::Identity() - 2*X23*X43 + X43*X43
    ) * Z42;

    return V_D + V_gauge;
}

row FZ34(
    row Z34, complex Z31, complex Z32, complex Z13, complex Z23, col Z41,
    col Z42, col Z43, row Z14, row Z24, matrix phi41, matrix phi42, matrix phi43,
    double X31, double X32, double X33, matrix X41, matrix X42, matrix X43)
{
    //V_D piece
    row V_D = complex(2,0) * (Z34*(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(Z13) - modsqr(Z23) - colabsqr(Z43) - c3/(g*g))
    - Z34*(Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34
    + commutator(phi41,phi41.adjoint()) + commutator(phi42,phi42.adjoint()) + commutator(phi43,phi43.adjoint()) - c4/(g*g) * matrix::Identity()));

    row V_gauge = complex(2,0) * Z34 * (
        X31*X31*matrix::Identity() - 2*X31*X41 + X41*X41 +
        X32*X32*matrix::Identity() - 2*X32*X42 + X42*X42 +
        X33*X33*matrix::Identity() - 2*X33*X43 + X43*X43
    );

    return V_D + V_gauge;
}

col FZ43(
    row Z34, complex Z31, complex Z32, complex Z13, complex Z23, col Z41,
    col Z42, col Z43, row Z14, row Z24, matrix phi41, matrix phi42, matrix phi43,
    double X31, double X32, double X33, matrix X41, matrix X42, matrix X43)
{
    //V_D piece
    col V_D = complex(2,0) * (- Z43*(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(Z13) - modsqr(Z23) - colabsqr(Z43) - c3/(g*g))
    + (Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34
    + commutator(phi41,phi41.adjoint()) + commutator(phi42,phi42.adjoint()) + commutator(phi43,phi43.adjoint()) - c4/(g*g) * matrix::Identity())*Z43);

    col V_gauge = complex(2,0) * (
        X31*X31*matrix::Identity() - 2*X31*X41 + X41*X41 +
        X32*X32*matrix::Identity() - 2*X32*X42 + X42*X42 +
        X33*X33*matrix::Identity() - 2*X33*X43 + X43*X43
    ) * Z43;

    return V_D + V_gauge;
}

// Now the X functions
// Scalars here 
double FX11(double X11, double X21, double X31, complex Z12, complex Z21, complex Z13, complex Z31, row Z14, col Z41, matrix X41)
{
    //V_gauge piece
    double V_gauge = 2 * ((X11 - X21)*(modsqr(Z12) + modsqr(Z21)) + (X11 - X31)*(modsqr(Z13) + modsqr(Z31)))
    + 2 * X11 * (rowabsqr(Z14) + colabsqr(Z41)) - 2 * (Z14 * X41 * Z14.adjoint()).trace().real() - 2 * (Z41.adjoint() * X41 * Z41).trace().real();
    return V_gauge;
}

double FX12(double X12, double X22, double X32, complex Z12, complex Z21, complex Z13, complex Z31, row Z14, col Z41, matrix X42)
{
    //V_gauge piece
    double V_gauge = 2 * ((X12 - X22)*(modsqr(Z12) + modsqr(Z21)) + (X12 - X32)*(modsqr(Z13) + modsqr(Z31)))
    + 2 * X12 * (rowabsqr(Z14) + colabsqr(Z41)) - 2 * (Z14 * X42 * Z14.adjoint()).trace().real() - 2 * (Z41.adjoint() * X42 * Z41).trace().real();
    return V_gauge;
}

double FX13(double X13, double X23, double X33, complex Z12, complex Z21, complex Z13, complex Z31, row Z14, col Z41, matrix X43)
{
    //V_gauge piece
    double V_gauge = 2 * ((X13 - X23)*(modsqr(Z12) + modsqr(Z21)) + (X13 - X33)*(modsqr(Z13) + modsqr(Z31)))
    + 2 * X13 * (rowabsqr(Z14) + colabsqr(Z41)) - 2 * (Z14 * X43 * Z14.adjoint()).trace().real() - 2 * (Z41.adjoint() * X43 * Z41).trace().real();
    return V_gauge;
}

double FX21(double X21, double X11, double X31, complex Z21, complex Z12, complex Z23, complex Z32, row Z24, col Z42, matrix X41)
{
    //V_gauge piece
    double V_gauge = 2 * ((X21 - X11)*(modsqr(Z21) + modsqr(Z12)) + (X21 - X31)*(modsqr(Z23) + modsqr(Z32)))
    + 2 * X21 * (rowabsqr(Z24) + colabsqr(Z42)) - 2 * (Z24 * X41 * Z24.adjoint()).trace().real() - 2 * (Z42.adjoint() * X41 * Z42).trace().real();
    return V_gauge; 
}

double FX22(double X22, double X12, double X32, complex Z21, complex Z12, complex Z23, complex Z32, row Z24, col Z42, matrix X42)
{
    double V_gauge = 2 * ((X22 - X12)*(modsqr(Z21) + modsqr(Z12)) + (X22 - X32)*(modsqr(Z23) + modsqr(Z32)))
    + 2 * X22 * (rowabsqr(Z24) + colabsqr(Z42)) - 2 * (Z24 * X42 * Z24.adjoint()).trace().real() - 2 * (Z42.adjoint() * X42 * Z42).trace().real();
    return V_gauge;
}

double FX23(double X23, double X13, double X33, complex Z21, complex Z12, complex Z23, complex Z32, row Z24, col Z42, matrix X43)
{
    double V_gauge = 2 * ((X23 - X13)*(modsqr(Z21) + modsqr(Z12)) + (X23 - X33)*(modsqr(Z23) + modsqr(Z32)))
    + 2 * X23 * (rowabsqr(Z24) + colabsqr(Z42)) - 2 * (Z24 * X43 * Z24.adjoint()).trace().real() - 2 * (Z42.adjoint() * X43 * Z42).trace().real();
    return V_gauge;
}

double FX31(double X31, double X11, double X21, complex Z31, complex Z13, complex Z32, complex Z23, row Z34, col Z43, matrix X41)
{
    double V_gauge = 2 * ((X31 - X11)*(modsqr(Z31) + modsqr(Z13)) + (X31 - X21)*(modsqr(Z32) + modsqr(Z23)))
    + 2 * X31 * (rowabsqr(Z34) + colabsqr(Z43)) - 2 * (Z34 * X41 * Z34.adjoint()).trace().real() - 2 * (Z43.adjoint() * X41 * Z43).trace().real();
    return V_gauge;
}

double FX32(double X32, double X12, double X22, complex Z31, complex Z13, complex Z32, complex Z23, row Z34, col Z43, matrix X42)
{
    double V_gauge = 2 * ((X32 - X12)*(modsqr(Z31) + modsqr(Z13)) + (X32 - X22)*(modsqr(Z32) + modsqr(Z23)))
    + 2 * X32 * (rowabsqr(Z34) + colabsqr(Z43)) - 2 * (Z34 * X42 * Z34.adjoint()).trace().real() - 2 * (Z43.adjoint() * X42 * Z43).trace().real();
    return V_gauge;
}

double FX33(double X33, double X13, double X23, complex Z31, complex Z13, complex Z32, complex Z23, row Z34, col Z43, matrix X43)
{
    double V_gauge = 2 * ((X33 - X13)*(modsqr(Z31) + modsqr(Z13)) + (X33 - X23)*(modsqr(Z32) + modsqr(Z23)))
    + 2 * X33 * (rowabsqr(Z34) + colabsqr(Z43)) - 2 * (Z34 * X43 * Z34.adjoint()).trace().real() - 2 * (Z43.adjoint() * X43 * Z43).trace().real();
    return V_gauge;
} 

/*
    If I compute Z4l*z4l.adjoint() and Zl4.adjoint()*Zl4 in the update function then it saves computation time in 
    the F's, just call them as arguments
*/

matrix FX41(
    double X11, double X21, double X31, 
    matrix X41, matrix X42, matrix X43, matrix phi41, matrix phi42, matrix phi43,
    matrix M41, matrix M42, matrix M43, matrix M14, matrix M24, matrix M34)
{
    // Working with V_gauge {splitting into 3 parts}
    // Store the matrices formed from Z4l*Z4l.adjoint() and Zl4.adjoint()Zl4 and calling them M4l and Ml4 respectively
    
    matrix V1_gauge = (
        anticommutator(M41,X41) - 2 * X11 *(M41 + M14) + anticommutator(M14,X41) +
        anticommutator(M42,X41) - 2 * X21 *(M42 + M24) + anticommutator(M24,X41) + 
        anticommutator(M43,X41) - 2 * X31 *(M43 + M34) + anticommutator(M34,X41)
    );

    // Storing some things computations for the second part calling them P's
    
    matrix P1 = commutator(phi41, commutator(phi41.adjoint(), X41)); 
    matrix P2 = commutator(phi42, commutator(phi42.adjoint(), X41));
    matrix P3 = commutator(phi43, commutator(phi43.adjoint(), X41));

    matrix V2_gauge = (
        (P1 + P1.adjoint()) + (P2 + P2.adjoint()) + (P3 + P3.adjoint())
    );

    // Now the third part
    matrix V3_gauge = (
        commutator(commutator(X41, X41), X41) + commutator(commutator(X41, X42), X42) + commutator(commutator(X41, X43), X43)
    );

    return V1_gauge + V2_gauge + V3_gauge;
}

matrix FX42(
    double X12, double X22, double X32, 
    matrix X41, matrix X42, matrix X43, matrix phi41, matrix phi42, matrix phi43,
    matrix M41, matrix M42, matrix M43, matrix M14, matrix M24, matrix M34
)
{
    // V1_gauge 
    matrix V1_gauge = (
        anticommutator(M41,X42) - 2 * X12 *(M41 + M14) + anticommutator(M14,X42) +
        anticommutator(M42,X42) - 2 * X22 *(M42 + M24) + anticommutator(M24,X42) +
        anticommutator(M43,X42) - 2 * X32 *(M43 + M34) + anticommutator(M34,X42)
    );

    // V2_gauge
    matrix P1 = commutator(phi41, commutator(phi41.adjoint(), X42)); 
    matrix P2 = commutator(phi42, commutator(phi42.adjoint(), X42));
    matrix P3 = commutator(phi43, commutator(phi43.adjoint(), X42));
    
    matrix V2_gauge = (
        (P1 + P1.adjoint()) + (P2 + P2.adjoint()) + (P3 + P3.adjoint())
    );
    // V3_gauge
    matrix V3_gauge = (
        commutator(commutator(X42,X41), X41) + commutator(commutator(X42, X42), X42) + commutator(commutator(X42, X43), X43)
    );

    return V1_gauge + V2_gauge + V3_gauge;
}

matrix FX43(
    double X13, double X23, double X33,
    matrix X41, matrix X42, matrix X43, matrix phi41, matrix phi42, matrix phi43,
    matrix M41, matrix M42, matrix M43, matrix M14, matrix M24, matrix M34
)
{
    // V1_gauge
    matrix V1_gauge = (
        anticommutator(M41,X43) - 2 * X13 *(M41 + M14) + anticommutator(M14,X43) +
        anticommutator(M42,X43) - 2 * X23 *(M42 + M24) + anticommutator(M24,X43) +
        anticommutator(M43,X43) - 2 * X33 *(M43 + M34) + anticommutator(M34,X43)
    );

    // V2_gauge
    matrix P1 = commutator(phi41, commutator(phi41.adjoint(), X43)); 
    matrix P2 = commutator(phi42, commutator(phi42.adjoint(), X43));
    matrix P3 = commutator(phi43, commutator(phi43.adjoint(), X43));

    matrix V2_gauge = (
        (P1 + P1.adjoint()) + (P2 + P2.adjoint()) + (P3 + P3.adjoint())
    );

    // V3_gauge
    matrix V3_gauge = (
        commutator(commutator(X43, X41), X41) + commutator(commutator(X43, X42), X42) + commutator(commutator(X43, X43), X43)
    );
    
    return V1_gauge + V2_gauge + V3_gauge;
}

matrix FX44(matrix X45, col Z41, col Z42, col Z43, row Z14, row Z24, row Z34, matrix phi41, matrix phi42, matrix phi43, double c4)
{
    //V_D part 
    matrix V_D = complex(0,2)*commutator(gamma(Z41, Z42, Z43, Z14, Z24, Z34, phi41, phi42, phi43, c4),X45);
    
    // V1_gauge
    matrix V1_gauge = (
        anticommutator(M41,X44) - 2 * X14 *(M41 + M14) + anticommutator(M14,X44) +
        anticommutator(M42,X44) - 2 * X24 *(M42 + M24) + anticommutator(M24,X44) +
        anticommutator(M43,X44) - 2 * X34 *(M43 + M34) + anticommutator(M34,X44)
    );
    
    // V2_gauge
    matrix P1 = commutator(X41, commutator(X41, phi41));
    matrix P2 = commutator(X42, commutator(X42, phi42));
    matrix P3 = commutator(X43, commutator(X43, phi43));

    matrix V2_gauge = (
        (P1 + P1.adjoint()) + (P2 + P2.adjoint()) + (P3 + P3.adjoint()) 
    ); 

    matrix V3_gauge = (
        commutator(commutator(X44, X41), X41) + commutator(commutator(X44,X42), X42) + commutator(commutator(X44,X43), X43)
    );
    return V_D + V1_gauge + V2_gauge;
}

matrix FX45(matrix X44, col Z41, col Z42, col Z43, row Z14, row Z24, row Z34, matrix phi41, matrix phi42, matrix phi43, double c4)
{
    //V_D part 
    matrix vd = complex(0,2)*commutator(X44, gamma(Z41, Z42, Z43, Z14, Z24, Z34, phi41, phi42, phi43, c4));
    return vd;
}

matrix FX46(matrix X47, col Z41, col Z42, col Z43, row Z14, row Z24, row Z34, matrix phi41, matrix phi42, matrix phi43, double c4)
{
    //V_D part
    matrix vd = complex(0,2)*commutator(gamma(Z41, Z42, Z43, Z14, Z24, Z34, phi41, phi42, phi43, c4),X47);
    return vd;
}

matrix FX47(matrix X46, col Z41, col Z42, col Z43, row Z14, row Z24, row Z34, matrix phi41, matrix phi42, matrix phi43, double c4)
{
    //V_D part 
    matrix vd = complex(0,2)*commutator(X46, gamma(Z41, Z42, Z43, Z14, Z24, Z34, phi41, phi42, phi43, c4));
    return vd;
}

matrix FX48(matrix X49, col Z41, col Z42, col Z43, row Z14, row Z24, row Z34, matrix phi41, matrix phi42, matrix phi43, double c4)
{
    //V_D part
    matrix vd = complex(0,2)*commutator(gamma(Z41, Z42, Z43, Z14, Z24, Z34, phi41, phi42, phi43, c4),X49);
    return vd;
}

matrix FX49(matrix X48, col Z41, col Z42, col Z43, row Z14, row Z24, row Z34, matrix phi41, matrix phi42, matrix phi43, double c4)
{
    //V_D part 
    matrix vd = complex(0,2)*commutator(X48, gamma(Z41, Z42, Z43, Z14, Z24, Z34, phi41, phi42, phi43, c4));
    return vd;
}

// Define the Update function... Huge scope for error here make sure that all of the paramaters of the force functions are in the 
// Correct order or else they will give sneaky errors {need to add x-coordinates and update the phi's}

void update(
    double dt, double* c4,
    double* V11, double* V12, double* V13, double* V14, double* V15, double* V16, double* V17, double* V18, double* V19,
    double* V21, double* V22, double* V23, double* V24, double* V25, double* V26, double* V27, double* V28, double* V29,
    double* V31, double* V32, double* V33, double* V34, double* V35, double* V36, double* V37, double* V38, double* V39,
    matrix* V41, matrix* V42, matrix* V43, matrix* V44, matrix* V45, matrix* V46, matrix* V47, matrix* V48, matrix* V49,
    double* X11, double* X12, double* X13, double* X14, double* X15, double* X16, double* X17, double* X18, double* X19,
    double* X21, double* X22, double* X23, double* X24, double* X25, double* X26, double* X27, double* X28, double* X29,
    double* X31, double* X32, double* X33, double* X34, double* X35, double* X36, double* X37, double* X38, double* X39,
    matrix* X41, matrix* X42, matrix* X43, matrix* X44, matrix* X45, matrix* X46, matrix* X47, matrix* X48, matrix* X49,
    double* X11_n, double* X12_n, double* X13_n, double* X14_n, double* X15_n, double* X16_n, double* X17_n, double* X18_n, double* X19_n,
    double* X21_n, double* X22_n, double* X23_n, double* X24_n, double* X25_n, double* X26_n, double* X27_n, double* X28_n, double* X29_n,
    double* X31_n, double* X32_n, double* X33_n, double* X34_n, double* X35_n, double* X36_n, double* X37_n, double* X38_n, double* X39_n,
    matrix* X41_n, matrix* X42_n, matrix* X43_n, matrix* X44_n, matrix* X45_n, matrix* X46_n, matrix* X47_n, matrix* X48_n, matrix* X49_n,
    matrix* phi41, matrix* phi42, matrix* phi43, matrix* phi41_n, matrix* phi42_n, matrix* phi43_n,
    complex* Z12, complex* Z12_n, complex* Z21, complex* Z21_n, complex* Z13, complex* Z13_n, complex* Z31, complex* Z31_n,
    complex* Z23, complex* Z23_n, complex* Z32, complex* Z32_n, col* Z41, col* Z41_n, col* Z42, col* Z42_n, col* Z43,
    col* Z43_n, row* Z14, row* Z14_n, row* Z24, row* Z24_n, row* Z34, row* Z34_n,
    complex* U12, complex* U13, complex* U21, complex* U23, complex* U31, complex* U32,
    row* U14, row* U24, row* U34, col* U41, col* U42, col* U43,
    matrix* M41, matrix* M42, matrix* M43, matrix* M14, matrix* M24, matrix* M34,
    matrix* M41_n, matrix* M42_n, matrix* M43_n, matrix* M14_n, matrix* M24_n, matrix* M34_n)
{
    
    // Updating the "Z_n coordinates"
    *Z12_n = *Z12 + *U12*dt - 0.5*FZ12(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z24, *Z32, *Z42, *X11, *X12, *X13, *X21, *X22, *X23)*pow(dt,2);
    *Z13_n = *Z13 + *U13*dt - 0.5*FZ13(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z34, *Z32, *Z43, *X11, *X12, *X13, *X31, *X32, *X33)*pow(dt,2);
    *Z21_n = *Z21 + *U21*dt - 0.5*FZ21(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z24, *Z32, *Z42, *X11, *X12, *X13, *X21, *X22, *X23)*pow(dt,2);
    *Z23_n = *Z23 + *U23*dt - 0.5*FZ23(*Z23, *Z21, *Z24, *Z12, *Z32, *Z42, *Z31, *Z34, *Z13, *Z43, *X21, *X22, *X23, *X31, *X32, *X33)*pow(dt,2);
    *Z31_n = *Z31 + *U31*dt - 0.5*FZ31(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z34, *Z32, *Z43, *X11, *X12, *X13, *X31, *X32, *X33)*pow(dt,2);
    *Z32_n = *Z32 + *U32*dt - 0.5*FZ32(*Z23, *Z21, *Z24, *Z12, *Z32, *Z42, *Z31, *Z34, *Z13, *Z43, *X21, *X22, *X23, *X31, *X32, *X33)*pow(dt,2);
    
    *Z14_n = *Z14 + *U14*dt - 0.5*FZ14(*Z14, *Z12, *Z13, *Z21, *Z31, *Z41, *Z42, *Z43, *Z24, *Z34, *phi41, *phi42, *phi43, *X11, *X12, *X13, *X41, *X42, *X43)*pow(dt,2);
    *Z24_n = *Z24 + *U24*dt - 0.5*FZ24(*Z24, *Z21, *Z23, *Z12, *Z32, *Z41, *Z42, *Z43, *Z14, *Z34, *phi41, *phi42, *phi43, *X21, *X22, *X23, *X41, *X42, *X43)*pow(dt,2);
    *Z34_n = *Z34 + *U34*dt - 0.5*FZ34(*Z34, *Z31, *Z32, *Z13, *Z23, *Z41, *Z42, *Z43, *Z14, *Z24, *phi41, *phi42, *phi43, *X31, *X32, *X33, *X41, *X42, *X43)*pow(dt,2);
    
    *Z41_n = *Z41 + *U41*dt - 0.5*FZ41(*Z14, *Z12, *Z13, *Z21, *Z31, *Z41, *Z42, *Z43, *Z24, *Z34, *phi41, *phi42, *phi43, *X11, *X12, *X13, *X41, *X42, *X43)*pow(dt,2);
    *Z42_n = *Z42 + *U42*dt - 0.5*FZ42(*Z24, *Z21, *Z23, *Z12, *Z32, *Z41, *Z42, *Z43, *Z14, *Z34, *phi41, *phi42, *phi43, *X21, *X22, *X23, *X41, *X42, *X43)*pow(dt,2);
    *Z43_n = *Z43 + *U43*dt - 0.5*FZ43(*Z34, *Z31, *Z32, *Z13, *Z23, *Z41, *Z42, *Z43, *Z14, *Z24, *phi41, *phi42, *phi43, *X31, *X32, *X33, *X41, *X42, *X43)*pow(dt,2);
    
    // Update these M coordinates 
    *M41_n = (*Z41_n) * (*Z41_n).adjoint(); 
    *M42_n = (*Z42_n) * (*Z42_n).adjoint(); 
    *M43_n = (*Z43_n) * (*Z43_n).adjoint();
    *M14_n = (*Z14_n).adjoint() * (*Z14_n); 
    *M24_n = (*Z24_n).adjoint() * (*Z24_n); 
    *M34_n = (*Z34_n).adjoint() * (*Z34_n);

    // Update the X_n coordinates and the phi_n's here
    *X11_n = *X11 + *V11*dt - 0.5*FX11(*X11, *X21, *X31, *Z12, *Z21, *Z13, *Z31, *Z14, *Z41, *X41)*pow(dt,2);
    *X12_n = *X12 + *V12*dt - 0.5*FX12(*X12, *X22, *X32, *Z12, *Z21, *Z13, *Z31, *Z14, *Z41, *X42)*pow(dt,2);
    *X13_n = *X13 + *V13*dt - 0.5*FX13(*X13, *X23, *X33, *Z12, *Z21, *Z13, *Z31, *Z14, *Z41, *X43)*pow(dt,2);
    *X14_n = *X14 + *V14*dt;
    *X15_n = *X15 + *V15*dt;
    *X16_n = *X16 + *V16*dt;
    *X17_n = *X17 + *V17*dt;
    *X18_n = *X18 + *V18*dt;
    *X19_n = *X19 + *V19*dt;
    *X21_n = *X21 + *V21*dt - 0.5*FX21(*X21, *X11, *X31, *Z21, *Z12, *Z23, *Z32, *Z24, *Z42, *X41)*pow(dt,2);
    *X22_n = *X22 + *V22*dt - 0.5*FX22(*X22, *X12, *X32, *Z21, *Z12, *Z23, *Z32, *Z24, *Z42, *X42)*pow(dt,2);
    *X23_n = *X23 + *V23*dt - 0.5*FX23(*X23, *X13, *X33, *Z21, *Z12, *Z23, *Z32, *Z24, *Z42, *X43)*pow(dt,2);
    *X24_n = *X24 + *V24*dt;
    *X25_n = *X25 + *V25*dt;
    *X26_n = *X26 + *V26*dt;
    *X27_n = *X27 + *V27*dt;
    *X28_n = *X28 + *V28*dt;
    *X29_n = *X29 + *V29*dt;

    *X31_n = *X31 + *V31*dt - 0.5*FX31(*X31, *X11, *X21, *Z31, *Z13, *Z32, *Z23, *Z34, *Z43, *X41)*pow(dt,2);
    *X32_n = *X32 + *V32*dt - 0.5*FX32(*X32, *X12, *X22, *Z31, *Z13, *Z32, *Z23, *Z34, *Z43, *X42)*pow(dt,2);
    *X33_n = *X33 + *V33*dt - 0.5*FX33(*X33, *X13, *X23, *Z31, *Z13, *Z32, *Z23, *Z34, *Z43, *X43)*pow(dt,2);
    *X34_n = *X34 + *V34*dt;
    *X35_n = *X35 + *V35*dt;
    *X36_n = *X36 + *V36*dt;
    *X37_n = *X37 + *V37*dt;
    *X38_n = *X38 + *V38*dt;
    *X39_n = *X39 + *V39*dt;

    *X41_n = *X41 + *V41*dt - 0.5*FX41(*X11, *X21, *X31, *X41, *X42, *X43, *phi41, *phi42, *phi43, *M41, *M42, *M43, *M14, *M24, *M34)*pow(dt,2);
    *X42_n = *X42 + *V42*dt - 0.5*FX42(*X12, *X22, *X32, *X41, *X42, *X43, *phi41, *phi42, *phi43, *M41, *M42, *M43, *M14, *M24, *M34)*pow(dt,2);
    *X43_n = *X43 + *V43*dt - 0.5*FX43(*X13, *X23, *X33, *X41, *X42, *X43, *phi41, *phi42, *phi43, *M41, *M42, *M43, *M14, *M24, *M34)*pow(dt,2);
    *X44_n = *X44 + *V44*dt - 0.5*FX44(*X45, *Z41, *Z42, *Z43, *Z14, *Z24, *Z34, *phi41, *phi42, *phi43, *c4)*pow(dt,2);
    *X45_n = *X45 + *V45*dt - 0.5*FX45(*X44, *Z41, *Z42, *Z43, *Z14, *Z24, *Z34, *phi41, *phi42, *phi43, *c4)*pow(dt,2);
    *X46_n = *X46 + *V46*dt - 0.5*FX46(*X47, *Z41, *Z42, *Z43, *Z14, *Z24, *Z34, *phi41, *phi42, *phi43, *c4)*pow(dt,2);
    *X47_n = *X47 + *V47*dt - 0.5*FX47(*X46, *Z41, *Z42, *Z43, *Z14, *Z24, *Z34, *phi41, *phi42, *phi43, *c4)*pow(dt,2);
    *X48_n = *X48 + *V48*dt - 0.5*FX48(*X49, *Z41, *Z42, *Z43, *Z14, *Z24, *Z34, *phi41, *phi42, *phi43, *c4)*pow(dt,2);
    *X49_n = *X49 + *V49*dt - 0.5*FX49(*X48, *Z41, *Z42, *Z43, *Z14, *Z24, *Z34, *phi41, *phi42, *phi43, *c4)*pow(dt,2);

    *phi41_n = *X44_n + complex(0,1) * *X45_n;
    *phi42_n = *X46_n + complex(0,1) * *X47_n;
    *phi43_n = *X48_n + complex(0,1) * *X49_n;

    // Update the U coordinates (The velocities of the Z's)
    *U12 = *U12 - 0.5*(FZ12(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z24, *Z32, *Z42, *X11, *X12, *X13, *X21, *X22, *X23) + FZ12(*Z12_n, *Z13_n, *Z14_n, *Z21_n, *Z31_n, *Z41_n, *Z23_n, *Z24_n, *Z32_n, *Z42_n, *X11_n, *X12_n, *X13_n, *X21_n, *X22_n, *X23_n))*dt;
    *U13 = *U13 - 0.5*(FZ13(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z34, *Z32, *Z43, *X11, *X12, *X13, *X31, *X32, *X33) + FZ13(*Z12_n, *Z13_n, *Z14_n, *Z21_n, *Z31_n, *Z41_n, *Z23_n, *Z34_n, *Z32_n, *Z43_n, *X11_n, *X12_n, *X13_n, *X31_n, *X32_n, *X33_n))*dt;
    *U21 = *U21 - 0.5*(FZ21(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z24, *Z32, *Z42, *X11, *X12, *X13, *X21, *X22, *X23) + FZ21(*Z12_n, *Z13_n, *Z14_n, *Z21_n, *Z31_n, *Z41_n, *Z23_n, *Z24_n, *Z32_n, *Z42_n, *X11_n, *X12_n, *X13_n, *X21_n, *X22_n, *X23_n))*dt;
    *U23 = *U23 - 0.5*(FZ23(*Z23, *Z21, *Z24, *Z12, *Z32, *Z42, *Z31, *Z34, *Z13, *Z43, *X21, *X22, *X23, *X31, *X32, *X33) + FZ23(*Z23_n, *Z21_n, *Z24_n, *Z12_n, *Z32_n, *Z42_n, *Z31_n, *Z34_n, *Z13_n, *Z43_n, *X21_n, *X22_n, *X23_n, *X31_n, *X32_n, *X33_n))*dt; 
    *U31 = *U31 - 0.5*(FZ31(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z34, *Z32, *Z43, *X11, *X12, *X13, *X31, *X32, *X33) + FZ31(*Z12_n, *Z13_n, *Z14_n, *Z21_n, *Z31_n, *Z41_n, *Z23_n, *Z34_n, *Z32_n, *Z43_n, *X11_n, *X12_n, *X13_n, *X31_n, *X32_n, *X33_n))*dt;
    *U32 = *U32 - 0.5*(FZ32(*Z23, *Z21, *Z24, *Z12, *Z32, *Z42, *Z31, *Z34, *Z13, *Z43, *X21, *X22, *X23, *X31, *X32, *X33) + FZ32(*Z23_n, *Z21_n, *Z24_n, *Z12_n, *Z32_n, *Z42_n, *Z31_n, *Z34_n, *Z13_n, *Z43_n, *X21_n, *X22_n, *X23_n, *X31_n, *X32_n, *X33_n))*dt;         
    *U14 = *U14 - 0.5*(FZ14(*Z14, *Z12, *Z13, *Z21, *Z31, *Z41, *Z42, *Z43, *Z24, *Z34, *phi41, *phi42, *phi43, *X11, *X12, *X13, *X41, *X42, *X43) + FZ14(*Z14_n, *Z12_n, *Z13_n, *Z21_n, *Z31_n, *Z41_n, *Z42_n, *Z43_n, *Z24_n, *Z34_n, *phi41_n, *phi42_n, *phi43_n, *X11_n, *X12_n, *X13_n, *X41_n, *X42_n, *X43_n))*dt;
    *U24 = *U24 - 0.5*(FZ24(*Z24, *Z21, *Z23, *Z12, *Z32, *Z41, *Z42, *Z43, *Z14, *Z34, *phi41, *phi42, *phi43, *X21, *X22, *X23, *X41, *X42, *X43) + FZ24(*Z24_n, *Z21_n, *Z23_n, *Z12_n, *Z32_n, *Z41_n, *Z42_n, *Z43_n, *Z14_n, *Z34_n, *phi41_n, *phi42_n, *phi43_n, *X21_n, *X22_n, *X23_n, *X41_n, *X42_n, *X43_n))*dt;
    *U34 = *U34 - 0.5*(FZ34(*Z34, *Z31, *Z32, *Z13, *Z23, *Z41, *Z42, *Z43, *Z14, *Z24, *phi41, *phi42, *phi43, *X31, *X32, *X33, *X41, *X42, *X43) + FZ34(*Z34_n, *Z31_n, *Z32_n, *Z13_n, *Z23_n, *Z41_n, *Z42_n, *Z43_n, *Z14_n, *Z24_n, *phi41_n, *phi42_n, *phi43_n, *X31_n, *X32_n, *X33_n, *X41_n, *X42_n, *X43_n))*dt;
    *U41 = *U41 - 0.5*(FZ41(*Z14, *Z12, *Z13, *Z21, *Z31, *Z41, *Z42, *Z43, *Z24, *Z34, *phi41, *phi42, *phi43, *X11, *X12, *X13, *X41, *X42, *X43) + FZ41(*Z14_n, *Z12_n, *Z13_n, *Z21_n, *Z31_n, *Z41_n, *Z42_n, *Z43_n, *Z24_n, *Z34_n, *phi41_n, *phi42_n, *phi43_n, *X11_n, *X12_n, *X13_n, *X41_n, *X42_n, *X43_n))*dt;
    *U42 = *U42 - 0.5*(FZ42(*Z24, *Z21, *Z23, *Z12, *Z32, *Z41, *Z42, *Z43, *Z14, *Z34, *phi41, *phi42, *phi43, *X21, *X22, *X23, *X41, *X42, *X43) + FZ42(*Z24_n, *Z21_n, *Z23_n, *Z12_n, *Z32_n, *Z41_n, *Z42_n, *Z43_n, *Z14_n, *Z34_n, *phi41_n, *phi42_n, *phi43_n, *X21_n, *X22_n, *X23_n, *X41_n, *X42_n, *X43_n))*dt;
    *U43 = *U43 - 0.5*(FZ43(*Z34, *Z31, *Z32, *Z13, *Z23, *Z41, *Z42, *Z43, *Z14, *Z24, *phi41, *phi42, *phi43, *X31, *X32, *X33, *X41, *X42, *X43) + FZ43(*Z34_n, *Z31_n, *Z32_n, *Z13_n, *Z23_n, *Z41_n, *Z42_n, *Z43_n, *Z14_n, *Z24_n, *phi41_n, *phi42_n, *phi43_n, *X31_n, *X32_n, *X33_n, *X41_n, *X42_n, *X43_n))*dt;

    // Update the V's here (The velocities of the X's)
    *V11 = *V11 - 0.5*(FX11(*X11, *X21, *X31, *Z12, *Z21, *Z13, *Z31, *Z14, *Z41, *X41) + FX11(*X11_n, *X21_n, *X31_n, *Z12_n, *Z21_n, *Z13_n, *Z31_n, *Z14_n, *Z41_n, *X41_n))*dt;
    *V12 = *V12 - 0.5*(FX12(*X12, *X22, *X32, *Z12, *Z21, *Z13, *Z31, *Z14, *Z41, *X42) + FX12(*X12_n, *X22_n, *X32_n, *Z12_n, *Z21_n, *Z13_n, *Z31_n, *Z14_n, *Z41_n, *X42_n))*dt; 
    *V13 = *V13 - 0.5*(FX13(*X13, *X23, *X33, *Z12, *Z21, *Z13, *Z31, *Z14, *Z41, *X43) + FX13(*X13_n, *X23_n, *X33_n, *Z12_n, *Z21_n, *Z13_n, *Z31_n, *Z14_n, *Z41_n, *X43_n))*dt; 
    *V14 = *V14; 
    *V15 = *V15; 
    *V16 = *V16; 
    *V17 = *V17; 
    *V18 = *V18; 
    *V19 = *V19;  
    *V21 = *V21 - 0.5*(FX21(*X21, *X11, *X31, *Z21, *Z12, *Z23, *Z32, *Z24, *Z42, *X41) + FX21(*X21_n, *X11_n, *X31_n, *Z21_n, *Z12_n, *Z23_n, *Z32_n, *Z24_n, *Z42_n, *X41_n))*dt;
    *V22 = *V22 - 0.5*(FX22(*X22, *X12, *X32, *Z21, *Z12, *Z23, *Z32, *Z24, *Z42, *X42) + FX22(*X22_n, *X12_n, *X32_n, *Z21_n, *Z12_n, *Z23_n, *Z32_n, *Z24_n, *Z42_n, *X42_n))*dt; 
    *V23 = *V23 - 0.5*(FX23(*X23, *X13, *X33, *Z21, *Z12, *Z23, *Z32, *Z24, *Z42, *X43) + FX23(*X23_n, *X13_n, *X33_n, *Z21_n, *Z12_n, *Z23_n, *Z32_n, *Z24_n, *Z42_n, *X43_n))*dt; 
    *V24 = *V24; 
    *V25 = *V25; 
    *V26 = *V26; 
    *V27 = *V27; 
    *V28 = *V28; 
    *V29 = *V29;

    *V31 = *V31 - 0.5*(FX31(*X31, *X11, *X21, *Z31, *Z13, *Z32, *Z23, *Z34, *Z43, *X41) + FX31(*X31_n, *X11_n, *X21_n, *Z31_n, *Z13_n, *Z32_n, *Z23_n, *Z34_n, *Z43_n, *X41_n))*dt;
    *V32 = *V32 - 0.5*(FX32(*X32, *X12, *X22, *Z31, *Z13, *Z32, *Z23, *Z34, *Z43, *X42) + FX32(*X32_n, *X12_n, *X22_n, *Z31_n, *Z13_n, *Z32_n, *Z23_n, *Z34_n, *Z43_n, *X42_n))*dt; 
    *V33 = *V33 - 0.5*(FX33(*X33, *X13, *X23, *Z31, *Z13, *Z32, *Z23, *Z34, *Z43, *X43) + FX33(*X33_n, *X13_n, *X23_n, *Z31_n, *Z13_n, *Z32_n, *Z23_n, *Z34_n, *Z43_n, *X43_n))*dt; 
    *V34 = *V34; 
    *V35 = *V35; 
    *V36 = *V36; 
    *V37 = *V37; 
    *V38 = *V38; 
    *V39 = *V39;

    *V41 = *V41 - 0.5*(FX41(*X11, *X21, *X31, *X41, *X42, *X43, *phi41, *phi42, *phi43, *M41, *M42, *M43, *M14, *M24, *M34) + FX41(*X11_n, *X21_n, *X31_n, *X41_n, *X42_n, *X43_n, *phi41_n, *phi42_n, *phi43_n, *M41_n, *M42_n, *M43_n, *M14_n, *M24_n, *M34_n))*dt;
    *V42 = *V42 - 0.5*(FX42(*X12, *X22, *X32, *X41, *X42, *X43, *phi41, *phi42, *phi43, *M41, *M42, *M43, *M14, *M24, *M34) + FX42(*X12_n, *X22_n, *X32_n, *X41_n, *X42_n, *X43_n, *phi41_n, *phi42_n, *phi43_n, *M41_n, *M42_n, *M43_n, *M14_n, *M24_n, *M34_n))*dt;
    *V43 = *V43 - 0.5*(FX43(*X13, *X23, *X33, *X41, *X42, *X43, *phi41, *phi42, *phi43, *M41, *M42, *M43, *M14, *M24, *M34) + FX43(*X13_n, *X23_n, *X33_n, *X41_n, *X42_n, *X43_n, *phi41_n, *phi42_n, *phi43_n, *M41_n, *M42_n, *M43_n, *M14_n, *M24_n, *M34_n))*dt;
    *V44 = *V44 - 0.5*(FX44(*X45, *Z41, *Z42, *Z43, *Z14, *Z24, *Z34, *phi41, *phi42, *phi43, *c4) + FX44(*X45_n, *Z41_n, *Z42_n, *Z43_n, *Z14_n, *Z24_n, *Z34_n, *phi41_n, *phi42_n, *phi43_n, *c4))*dt;
    *V45 = *V45 - 0.5*(FX45(*X44, *Z41, *Z42, *Z43, *Z14, *Z24, *Z34, *phi41, *phi42, *phi43, *c4) + FX45(*X44_n, *Z41_n, *Z42_n, *Z43_n, *Z14_n, *Z24_n, *Z34_n, *phi41_n, *phi42_n, *phi43_n, *c4))*dt;
    *V46 = *V46 - 0.5*(FX46(*X47, *Z41, *Z42, *Z43, *Z14, *Z24, *Z34, *phi41, *phi42, *phi43, *c4) + FX46(*X47_n, *Z41_n, *Z42_n, *Z43_n, *Z14_n, *Z24_n, *Z34_n, *phi41_n, *phi42_n, *phi43_n, *c4))*dt;
    *V47 = *V47 - 0.5*(FX47(*X46, *Z41, *Z42, *Z43, *Z14, *Z24, *Z34, *phi41, *phi42, *phi43, *c4) + FX47(*X46_n, *Z41_n, *Z42_n, *Z43_n, *Z14_n, *Z24_n, *Z34_n, *phi41_n, *phi42_n, *phi43_n, *c4))*dt;
    *V48 = *V48 - 0.5*(FX48(*X49, *Z41, *Z42, *Z43, *Z14, *Z24, *Z34, *phi41, *phi42, *phi43, *c4) + FX48(*X49_n, *Z41_n, *Z42_n, *Z43_n, *Z14_n, *Z24_n, *Z34_n, *phi41_n, *phi42_n, *phi43_n, *c4))*dt;
    *V49 = *V49 - 0.5*(FX49(*X48, *Z41, *Z42, *Z43, *Z14, *Z24, *Z34, *phi41, *phi42, *phi43, *c4) + FX49(*X48_n, *Z41_n, *Z42_n, *Z43_n, *Z14_n, *Z24_n, *Z34_n, *phi41_n, *phi42_n, *phi43_n, *c4))*dt;

    // The final updation
    // Z's go first
    *Z12 = *Z12_n;
    *Z13 = *Z13_n;
    *Z21 = *Z21_n;
    *Z23 = *Z23_n;
    *Z31 = *Z31_n;
    *Z32 = *Z32_n;
    *Z14 = *Z14_n;
    *Z24 = *Z24_n;
    *Z34 = *Z34_n;
    *Z41 = *Z41_n;
    *Z42 = *Z42_n;
    *Z43 = *Z43_n;

    // Then the X's
    *X11 = *X11_n;
    *X12 = *X12_n;
    *X13 = *X13_n;
    *X14 = *X14_n;
    *X15 = *X15_n;
    *X16 = *X16_n;
    *X17 = *X17_n;
    *X18 = *X18_n;
    *X19 = *X19_n;

    *X21 = *X21_n;
    *X22 = *X22_n;
    *X23 = *X23_n;
    *X24 = *X24_n;
    *X25 = *X25_n;
    *X26 = *X26_n;
    *X27 = *X27_n;
    *X28 = *X28_n;
    *X29 = *X29_n;

    *X31 = *X31_n;
    *X32 = *X32_n;
    *X33 = *X33_n;
    *X34 = *X34_n;
    *X35 = *X35_n;
    *X36 = *X36_n;
    *X37 = *X37_n;
    *X38 = *X38_n;
    *X39 = *X39_n;
    
    *X41 = *X41_n;
    *X42 = *X42_n;
    *X43 = *X43_n;
    *X44 = *X44_n;
    *X45 = *X45_n;
    *X46 = *X46_n;
    *X47 = *X47_n;
    *X48 = *X48_n;
    *X49 = *X49_n;

    // Finally the phi's
    *phi41 = *phi41_n;
    *phi42 = *phi42_n;
    *phi43 = *phi43_n;

    *M41 = *M41_n;
    *M42 = *M42_n;
    *M43 = *M43_n;
    *M14 = *M14_n;
    *M24 = *M24_n;
    *M34 = *M34_n;
}

complex Z12_sol[iterations/1000];
complex Z13_sol[iterations/1000];
complex Z21_sol[iterations/1000];
complex Z23_sol[iterations/1000];
complex Z31_sol[iterations/1000];
complex Z32_sol[iterations/1000];
row Z14_sol[iterations/1000];
row Z24_sol[iterations/1000];
row Z34_sol[iterations/1000];
col Z41_sol[iterations/1000];
col Z42_sol[iterations/1000];
col Z43_sol[iterations/1000];

int main()
{
    // Build the matrices such that ther are hermitian (The have already been declared but c++ only seems to like "for" loops when they are in main)
    for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                X41(i, j) = complex(gauss_dist(rng),gauss_dist(rng)), X41(j, i) = std:: conj(X41(i,j));
                X42(i, j) = complex(gauss_dist(rng),gauss_dist(rng)), X42(j, i) = std:: conj(X42(i,j));
                X43(i, j) = complex(gauss_dist(rng),gauss_dist(rng)), X43(j, i) = std:: conj(X43(i,j));
                X44(i, j) = complex(gauss_dist(rng),gauss_dist(rng)), X44(j, i) = std:: conj(X44(i,j));
                X45(i, j) = complex(gauss_dist(rng),gauss_dist(rng)), X45(j, i) = std:: conj(X45(i,j));
                X46(i, j) = complex(gauss_dist(rng),gauss_dist(rng)), X46(j, i) = std:: conj(X46(i,j));
                X47(i, j) = complex(gauss_dist(rng),gauss_dist(rng)), X47(j, i) = std:: conj(X47(i,j));
                X48(i, j) = complex(gauss_dist(rng),gauss_dist(rng)), X48(j, i) = std:: conj(X48(i,j));
                X49(i, j) = complex(gauss_dist(rng),gauss_dist(rng)), X49(j, i) = std:: conj(X49(i,j));
                if (i == j && i != N-1)
                {
                    X41(i,j) = complex(gauss_dist(rng),0);
                    X42(i,j) = complex(gauss_dist(rng),0);
                    X43(i,j) = complex(gauss_dist(rng),0);
                    X44(i,j) = complex(gauss_dist(rng),0);
                    X45(i,j) = complex(gauss_dist(rng),0);
                    X46(i,j) = complex(gauss_dist(rng),0);
                    X47(i,j) = complex(gauss_dist(rng),0);
                    X48(i,j) = complex(gauss_dist(rng),0);
                    X49(i,j) = complex(gauss_dist(rng),0);
                }
                if (i == N-1 && i == j)
                {
                    complex sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0, sum8 = 0, sum9 = 0;
                    for (int k = 0; k < N-1; k++)
                    {
                        sum1 += X41(k, k), sum2 += X42(k, k), sum3 += X43(k, k), sum4 += X44(k, k), sum5 += X45(k, k),
                        sum6 += X46(k, k), sum7 += X47(k, k), sum8 += X48(k, k), sum9 += X49(k, k);
                    }
                    X41(i, j) = -sum1, X42(i, j) = -sum2, X43(i, j) = -sum3, X44(i, j) = -sum4, X45(i, j) = -sum5,
                    X46(i, j) = -sum6, X47(i, j) = -sum7, X48(i, j) = -sum8, X49(i, j) = -sum9;
                }
            }
        }
    
    //Settingmatrices back to zero so that I can add them one at a time 
    
    // X41 = matrix::Zero();
    // X42 = matrix::Zero();
    // X43 = matrix::Zero();
    X44 = matrix::Zero();
    X45 = matrix::Zero();
    X46 = matrix::Zero();
    X47 = matrix::Zero();
    X48 = matrix::Zero();
    X49 = matrix::Zero();
    

    matrix phi41 = X44 + complex(0,1)*X45;
    matrix phi42 = X46 + complex(0,1)*X47;
    matrix phi43 = X48 + complex(0,1)*X49;
    matrix phi41_n = matrix::Zero();    
    matrix phi42_n = matrix::Zero();
    matrix phi43_n = matrix::Zero();

    for (int i = 0; i < iterations; i++)
    { 
        if (i % 1000 == 0)
        {
            Z12_sol[i/1000] = Z12;
            Z13_sol[i/1000] = Z13;
            Z21_sol[i/1000] = Z21;
            Z23_sol[i/1000] = Z23;
            Z31_sol[i/1000] = Z31;
            Z32_sol[i/1000] = Z32;
            Z14_sol[i/1000] = Z14;
            Z24_sol[i/1000] = Z24;
            Z34_sol[i/1000] = Z34;
            Z41_sol[i/1000] = Z41;
            Z42_sol[i/1000] = Z42;
            Z43_sol[i/1000] = Z43;
            
            std::cout << std::setprecision(10) <<K(
                U12, U13, U21, U23, U31, U32, U14, U24, U34, U41, U42, U43,
                V11, V12, V13, V14, V15, V16, V17, V18, V19,
                V21, V22, V23, V24, V25, V26, V27, V28, V29,
                V31, V32, V33, V34, V35, V36, V37, V38, V39,
                V41, V42, V43, V44, V45, V46, V47, V48, V49
                ) + 
                V_D(Z12, Z13, Z21, Z23, Z31, Z32, Z14, Z24, Z34, Z41, Z42, Z43, c1, c2, c3, c4, phi41, phi42, phi43) +
                V_gauge(Z12, Z13, Z21, Z23, Z31, Z32, Z14, Z24, Z34, Z41, Z42, Z43, X11, X12, X13, X21, X22, X23, X31, X32, X33,
                        X41, X42, X43, phi11, phi12, phi13, phi21, phi22, phi23, phi31, phi32, phi33, phi41, phi42, phi43) << std::endl;
            std::cout << "Time: " << (i*dt) << std::endl;   
        }

        update(
        dt, &c4,
        &V11, &V12, &V13, &V14, &V15, &V16, &V17, &V18, &V19,
        &V21, &V22, &V23, &V24, &V25, &V26, &V27, &V28, &V29,
        &V31, &V32, &V33, &V34, &V35, &V36, &V37, &V38, &V39,
        &V41, &V42, &V43, &V44, &V45, &V46, &V47, &V48, &V49,
        &X11, &X12, &X13, &X14, &X15, &X16, &X17, &X18, &X19,
        &X21, &X22, &X23, &X24, &X25, &X26, &X27, &X28, &X29,
        &X31, &X32, &X33, &X34, &X35, &X36, &X37, &X38, &X39,
        &X41, &X42, &X43, &X44, &X45, &X46, &X47, &X48, &X49,
        &X11_n, &X12_n, &X13_n, &X14_n, &X15_n, &X16_n, &X17_n, &X18_n, &X19_n,
        &X21_n, &X22_n, &X23_n, &X24_n, &X25_n, &X26_n, &X27_n, &X28_n, &X29_n,
        &X31_n, &X32_n, &X33_n, &X34_n, &X35_n, &X36_n, &X37_n, &X38_n, &X39_n,
        &X41_n, &X42_n, &X43_n, &X44_n, &X45_n, &X46_n, &X47_n, &X48_n, &X49_n,
        &phi41, &phi42, &phi43, &phi41_n, &phi42_n, &phi43_n,
        &Z12, &Z12_n, &Z21, &Z21_n, &Z13, &Z13_n, &Z31, &Z31_n,
        &Z23, &Z23_n, &Z32, &Z32_n, &Z41, &Z41_n, &Z42, &Z42_n, &Z43,   
        &Z43_n, &Z14, &Z14_n, &Z24, &Z24_n, &Z34, &Z34_n,
        &U12, &U13, &U21, &U23, &U31, &U32,
        &U14, &U24, &U34, &U41, &U42, &U43,
        &M41, &M42, &M43, &M14, &M24, &M34, &M41_n, &M42_n, &M43_n, &M14_n, &M24_n, &M34_n);    
    }

    return 0;
}