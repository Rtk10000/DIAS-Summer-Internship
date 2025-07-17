#include <iostream>
#include <random>
#include <ctime>
#include <string>
#include <complex>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "eigen/Eigen/Dense"

long double start = std::time(nullptr);

const long double g = 1;
const int N = 3;
const long double dt = 1e-4;
const int iterations = 1e5;

typedef std::complex<long double> complex;
typedef Eigen:: Matrix<complex, N, N> matrix;
typedef Eigen:: Matrix<complex, 1, N> row;
typedef Eigen:: Matrix<complex, N, 1> col;

long double colabsqr(col a)
{
    complex b = (a.adjoint()*a);
    return b.real();
}

long double rowabsqr(row a)
{
    complex b = (a*a.adjoint());
    return b.real();
}

long double modsqr(complex a)
{
    return pow(abs(a),2);
}

// Kinetic energy with the Z coordinates only 

long double K(
    complex U12, complex U13, complex U21, complex U23, complex U31, complex U32, row U14, row U24, row U34, col U41, col U42, col U43)
{
    long double scalar_sum = modsqr(U12) + modsqr(U13) + modsqr(U21) + modsqr(U23) + modsqr(U31) + modsqr(U32);
    long double row_sum = (U14.adjoint()*U14 + U24.adjoint()*U24 + U34.adjoint()*U34).trace().real();
    long double col_sum = (U41.adjoint()*U41 + U42.adjoint()*U42 + U43.adjoint()*U43).trace().real();
    return 0.5 * (scalar_sum + row_sum + col_sum);
}

long double V_D(
    complex Z12, complex Z13, complex Z21, complex Z23, complex Z31, complex Z32, row Z14, row Z24, row Z34, col Z41, col Z42, col Z43,
    long double c1, long double c2, long double c3, long double c4)
{
    // For k = 1
    long double sum1 = pow(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) - c1/g*g, 2);
    
    // For k = 2
    long double sum2 = pow(modsqr(Z21) + modsqr(Z23) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/g*g, 2);
    
    // For k = 3
    long double sum3 = pow(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(Z13) - modsqr(Z23) - colabsqr(Z43) - c3/g*g, 2);

    // For k = 4
    long double sum4 = ((
        Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34 - 
        (c4/g*g) * matrix::Identity()
        )*(
        Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34 - 
        (c4/g*g) * matrix::Identity()
        )).trace().real();    
        
    return 0.5 * (sum1 + sum2 + sum3 + sum4);
}
std:: mt19937 rng(std::time(nullptr));
std:: normal_distribution<long double> gauss_dist(0, 1);

// Initialise the Z variables only and let U be the \dot{Z} coordinates
/*
complex Z12 = complex(gauss_dist(rng),gauss_dist(rng));
complex Z13 = complex(gauss_dist(rng),gauss_dist(rng));
complex Z21 = complex(gauss_dist(rng),gauss_dist(rng));
complex Z23 = complex(gauss_dist(rng),gauss_dist(rng)); 
complex Z31 = complex(gauss_dist(rng),gauss_dist(rng)); 
complex Z32 = complex(gauss_dist(rng),gauss_dist(rng)); 

col Z41 = col::Random();
col Z42 = col::Random();
col Z43 = col::Random();
row Z14 = row::Random();
row Z24 = row::Random(); 
row Z34 = row::Random();
*/

complex Z12 = complex(0.207606,0.366979);
complex Z13 = complex(0.0365052,-0.344526);

complex Z21 = complex(0.0786534,-1.71609);
complex Z23 = complex(0.508537,1.50702);

complex Z31 = complex(0.547891,0.496717);
complex Z32 = complex(-0.937337,-1.09565);
//complex Z13 = complex(0,0);
//complex Z21 = complex(0,0);
//complex Z23 = complex(0,0);
//complex Z31 = complex(0,0);
//complex Z32 = complex(0,0);


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

// Initialise U
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

// Initialise c's as complex but they may be real ask Swapno if things don't work
/*
double c1 = gauss_dist(rng);
double c2 = gauss_dist(rng);
double c3 = gauss_dist(rng);
double c4 = gauss_dist(rng);
*/
double c1 = 2.5;
double c2 = 3.3;
double c3 = 1.2;
double c4 = -c1-c2-c3;

// Define the Energy 
double H(double K, double V_D)
{
    double E = K + V_D;
    return E;
}

// Define the Force functions (defined as in dropbox -\ddot(Zij) = FZij)

complex FZ12(
    complex Z12, complex Z13, row Z14, complex Z21, complex Z31, col Z41,
    complex Z23, row Z24, complex Z32, col Z42)
{
    return complex(2,0) * (Z12*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) - c1/(g*g))
    - Z12*(modsqr(Z21) + modsqr(Z23) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/(g*g)));
}

complex FZ21(
    complex Z12, complex Z13, row Z14, complex Z21, complex Z31, col Z41,
    complex Z23, row Z24, complex Z32, col Z42)
{
    return complex(2,0) * (- Z21*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) - c1/(g*g))
    + Z21*(modsqr(Z21) + modsqr(Z23) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/(g*g)));
}

complex FZ13(
    complex Z12, complex Z13, row Z14, complex Z21, complex Z31, col Z41,
    complex Z23, row Z34, complex Z32, col Z43)
{
    return complex(2,0) * (Z13*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) - c1/(g*g))
    - Z13*(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(Z13) - modsqr(Z23) - colabsqr(Z43) - c3/(g*g))); 
}

complex FZ31(
    complex Z12, complex Z13, row Z14, complex Z21, complex Z31, col Z41,
    complex Z23, row Z34, complex Z32, col Z43)
{
    return complex(2,0) * (- Z31*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) - c1/(g*g))
    + Z31*(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(Z13) - modsqr(Z23) - colabsqr(Z43) - c3/(g*g)));
}

complex FZ23(
    complex Z23, complex Z21, row Z24, complex Z12, complex Z32, col Z42,
    complex Z31, row Z34, complex Z13, col Z43)
{
    return complex(2,0) * (Z23*(modsqr(Z23) + modsqr(Z21) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/(g*g))
    - Z23*(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(Z13) - modsqr(Z23) - colabsqr(Z43) - c3/(g*g)));
}

complex FZ32(
    complex Z23, complex Z21, row Z24, complex Z12, complex Z32, col Z42,
    complex Z31, row Z34, complex Z13, col Z43)
{
    return complex(2,0) * (- Z32*(modsqr(Z23) + modsqr(Z21) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/(g*g))
    + Z32*(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(Z13) - modsqr(Z23) - colabsqr(Z43) - c3/(g*g)));
}

row FZ14(
    row Z14, complex Z12, complex Z13, complex Z21, complex Z31, col Z41, 
    col Z42, col Z43, row Z24, row Z34)
{
    return complex(2,0) * (Z14*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) -c1/(g*g))
    - Z14*(Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34
    - c4/(g*g) * matrix::Identity())); 
}

col FZ41(
    row Z14, complex Z12, complex Z13, complex Z21, complex Z31, col Z41, 
    col Z42, col Z43, row Z24, row Z34)
{
    return complex(2,0) * (- Z41*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) -c1/(g*g))
    + (Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34
    - c4/(g*g) * matrix::Identity())*Z41);
}

row FZ24(
    row Z24, complex Z21, complex Z23, complex Z12, complex Z32, col Z41,
    col Z42, col Z43, row Z14, row Z34)
{
    return complex(2,0) * (Z24*(modsqr(Z21) + modsqr(Z23) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/(g*g))
    - Z24*(Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34
    - c4/(g*g) * matrix::Identity()));
}

col FZ42(
    row Z24, complex Z21, complex Z23, complex Z12, complex Z32, col Z41,
    col Z42, col Z43, row Z14, row Z34)
{
    return complex(2,0) * (- Z42*(modsqr(Z21) + modsqr(Z23) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/(g*g))
    + (Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34
    - c4/(g*g) * matrix::Identity())*Z42);
}

row FZ34(
    row Z34, complex Z31, complex Z32, complex Z13, complex Z23, col Z41,
    col Z42, col Z43, row Z14, row Z24)
{
    return complex(2,0) * (Z34*(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(Z13) - modsqr(Z23) - colabsqr(Z43) - c3/(g*g))
    - Z34*(Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34
    - c4/(g*g) * matrix::Identity()));
}

col FZ43(
    row Z34, complex Z31, complex Z32, complex Z13, complex Z23, col Z41,
    col Z42, col Z43, row Z14, row Z24)
{
    return complex(2,0) * (- Z43*(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(Z13) - modsqr(Z23) - colabsqr(Z43) - c3/(g*g))
    + (Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34
    - c4/(g*g) * matrix::Identity())*Z43);
}

// Define the Update function... Huge scope for error here make sure that all of the paramaters of the force functions are in the 
// Correct order or else they will give unseen errors

void update(
    long double dt,
    complex* Z12, complex* Z12_n, complex* Z21, complex* Z21_n, complex* Z13, complex* Z13_n, complex* Z31, complex* Z31_n,
    complex* Z23, complex* Z23_n, complex* Z32, complex* Z32_n, col* Z41, col* Z41_n, col* Z42, col* Z42_n, col* Z43,
    col* Z43_n, row* Z14, row* Z14_n, row* Z24, row* Z24_n, row* Z34, row* Z34_n,
    complex* U12, complex* U13, complex* U21, complex* U23, complex* U31, complex* U32,
    row* U14, row* U24, row* U34, col* U41, col* U42, col* U43)
{
    *Z12_n = *Z12 + *U12*dt - complex(0.5, 0)*FZ12(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z24, *Z32, *Z42)*dt*dt;
    *Z13_n = *Z13 + *U13*dt - complex(0.5, 0)*FZ13(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z34, *Z32, *Z43)*dt*dt;
    *Z21_n = *Z21 + *U21*dt - complex(0.5, 0)*FZ21(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z24, *Z32, *Z42)*dt*dt;
    *Z23_n = *Z23 + *U23*dt - complex(0.5, 0)*FZ23(*Z23, *Z21, *Z24, *Z12, *Z32, *Z42, *Z31, *Z34, *Z13, *Z43)*dt*dt;
    *Z31_n = *Z31 + *U31*dt - complex(0.5, 0)*FZ31(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z34, *Z32, *Z43)*dt*dt;
    *Z32_n = *Z32 + *U32*dt - complex(0.5, 0)*FZ32(*Z23, *Z21, *Z24, *Z12, *Z32, *Z42, *Z31, *Z34, *Z13, *Z43)*dt*dt;
    *Z14_n = *Z14 + *U14*dt - complex(0.5, 0)*FZ14(*Z14, *Z12, *Z13, *Z21, *Z31, *Z41, *Z42, *Z43, *Z24, *Z34)*dt*dt;
    *Z24_n = *Z24 + *U24*dt - complex(0.5, 0)*FZ24(*Z24, *Z21, *Z23, *Z12, *Z32, *Z41, *Z42, *Z43, *Z14, *Z34)*dt*dt;
    *Z34_n = *Z34 + *U34*dt - complex(0.5, 0)*FZ34(*Z34, *Z31, *Z32, *Z13, *Z23, *Z41, *Z42, *Z43, *Z14, *Z24)*dt*dt;
    *Z41_n = *Z41 + *U41*dt - complex(0.5, 0)*FZ41(*Z14, *Z12, *Z13, *Z21, *Z31, *Z41, *Z42, *Z43, *Z24, *Z34)*dt*dt;
    *Z42_n = *Z42 + *U42*dt - complex(0.5, 0)*FZ42(*Z24, *Z21, *Z23, *Z12, *Z32, *Z41, *Z42, *Z43, *Z14, *Z34)*dt*dt;
    *Z43_n = *Z43 + *U43*dt - complex(0.5, 0)*FZ43(*Z34, *Z31, *Z32, *Z13, *Z23, *Z41, *Z42, *Z43, *Z14, *Z24)*dt*dt;

    *U12 = *U12 - complex(0.5, 0)*(FZ12(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z24, *Z32, *Z42) + FZ12(*Z12_n, *Z13_n, *Z14_n, *Z21_n, *Z31_n, *Z41_n, *Z23_n, *Z24_n, *Z32_n, *Z42_n))*dt;
    *U13 = *U13 - complex(0.5, 0)*(FZ13(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z34, *Z32, *Z43) + FZ13(*Z12_n, *Z13_n, *Z14_n, *Z21_n, *Z31_n, *Z41_n, *Z23_n, *Z34_n, *Z32_n, *Z43_n))*dt;
    *U21 = *U21 - complex(0.5, 0)*(FZ21(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z24, *Z32, *Z42) + FZ21(*Z12_n, *Z13_n, *Z14_n, *Z21_n, *Z31_n, *Z41_n, *Z23_n, *Z24_n, *Z32_n, *Z42_n))*dt;
    *U23 = *U23 - complex(0.5, 0)*(FZ23(*Z23, *Z21, *Z24, *Z12, *Z32, *Z42, *Z31, *Z34, *Z13, *Z43) + FZ23(*Z23_n, *Z21_n, *Z24_n, *Z12_n, *Z32_n, *Z42_n, *Z31_n, *Z34_n, *Z13_n, *Z43_n))*dt; 
    *U31 = *U31 - complex(0.5, 0)*(FZ31(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z34, *Z32, *Z43) + FZ31(*Z12_n, *Z13_n, *Z14_n, *Z21_n, *Z31_n, *Z41_n, *Z23_n, *Z34_n, *Z32_n, *Z43_n))*dt;
    *U32 = *U32 - complex(0.5, 0)*(FZ32(*Z23, *Z21, *Z24, *Z12, *Z32, *Z42, *Z31, *Z34, *Z13, *Z43) + FZ32(*Z23_n, *Z21_n, *Z24_n, *Z12_n, *Z32_n, *Z42_n, *Z31_n, *Z34_n, *Z13_n, *Z43_n))*dt;         
    *U14 = *U14 - complex(0.5, 0)*(FZ14(*Z14, *Z12, *Z13, *Z21, *Z31, *Z41, *Z42, *Z43, *Z24, *Z34) + FZ14(*Z14_n, *Z12_n, *Z13_n, *Z21_n, *Z31_n, *Z41_n, *Z42_n, *Z43_n, *Z24_n, *Z34_n))*dt;
    *U24 = *U24 - complex(0.5, 0)*(FZ24(*Z24, *Z21, *Z23, *Z12, *Z32, *Z41, *Z42, *Z43, *Z14, *Z34) + FZ24(*Z24_n, *Z21_n, *Z23_n, *Z12_n, *Z32_n, *Z41_n, *Z42_n, *Z43_n, *Z14_n, *Z34_n))*dt;
    *U34 = *U34 - complex(0.5, 0)*(FZ34(*Z34, *Z31, *Z32, *Z13, *Z23, *Z41, *Z42, *Z43, *Z14, *Z24) + FZ34(*Z34_n, *Z31_n, *Z32_n, *Z13_n, *Z23_n, *Z41_n, *Z42_n, *Z43_n, *Z14_n, *Z24_n))*dt;
    *U41 = *U41 - complex(0.5, 0)*(FZ41(*Z14, *Z12, *Z13, *Z21, *Z31, *Z41, *Z42, *Z43, *Z24, *Z34) + FZ41(*Z14_n, *Z12_n, *Z13_n, *Z21_n, *Z31_n, *Z41_n, *Z42_n, *Z43_n, *Z24_n, *Z34_n))*dt;
    *U42 = *U42 - complex(0.5, 0)*(FZ42(*Z24, *Z21, *Z23, *Z12, *Z32, *Z41, *Z42, *Z43, *Z14, *Z34) + FZ42(*Z24_n, *Z21_n, *Z23_n, *Z12_n, *Z32_n, *Z41_n, *Z42_n, *Z43_n, *Z14_n, *Z34_n))*dt;
    *U43 = *U43 - complex(0.5, 0)*(FZ43(*Z34, *Z31, *Z32, *Z13, *Z23, *Z41, *Z42, *Z43, *Z14, *Z24) + FZ43(*Z34_n, *Z31_n, *Z32_n, *Z13_n, *Z23_n, *Z41_n, *Z42_n, *Z43_n, *Z14_n, *Z24_n))*dt;

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

 
    row Z14;
    Z14(0,0) = complex(1,1);
    Z14(0,1) = complex(2,1);
    Z14(0,2) = complex(3,1);
    row Z24;
    Z24(0,0) = complex(1,2);
    Z24(0,1) = complex(2,2);
    Z24(0,2) = complex(3,2);
    row Z34;
    Z34(0,0) = complex(1,3);
    Z34(0,1) = complex(2,3);
    Z34(0,2) = complex(3,3);

    col Z41;
    Z41(0,0) = complex(1,1);
    Z41(1,0) = complex(1,2);
    Z41(2,0) = complex(1,3);
    col Z42;
    Z42(0,0) = complex(2,1);
    Z42(1,0) = complex(2,2);
    Z42(2,0) = complex(2,3);
    col Z43;
    Z43(0,0) = complex(3,1);
    Z43(1,0) = complex(3,2);
    Z43(2,0) = complex(3,3);
/*
    row Z14;
    Z14(0,0) = complex(0,0);
    Z14(0,1) = complex(0,0);
    Z14(0,2) = complex(0,0);
    row Z24;
    Z24(0,0) = complex(0,0);
    Z24(0,1) = complex(0,0);
    Z24(0,2) = complex(0,0);
    row Z34;
    Z34(0,0) = complex(0,0);
    Z34(0,1) = complex(0,0);
    Z34(0,2) = complex(0,0);

    col Z41;
    Z41(0,0) = complex(0,0);
    Z41(1,0) = complex(0,0);
    Z41(2,0) = complex(0,0);
    col Z42;
    Z42(0,0) = complex(0,0);
    Z42(1,0) = complex(0,0);
    Z42(2,0) = complex(0,0);
    col Z43;
    Z43(0,0) = complex(0,0);
    Z43(1,0) = complex(0,0);
    Z43(2,0) = complex(0,0);
*/

    std::cout << std::setprecision(15) << K(U12,U13,U21,U23,U31,U32,U14,U24,U34,U41,U42,U43)+V_D(Z12,Z13,Z21,Z23,Z31,Z32,Z14,Z24,Z34,Z41,Z42,Z43,c1,c2,c3,c4) << std::endl;
    std::cout << std::setprecision(15) << "Z12 = " << Z12;


    for (int i = 0; i < iterations; i++)
    {
        if (i % 1000 == 0)
        {
            /*
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
            */
            std::cout << '\n' << "H at" << i*dt << "seconds: ";
            std::cout << std::setprecision(15) << K(U12,U13,U21,U23,U31,U32,U14,U24,U34,U41,U42,U43)+V_D(Z12,Z13,Z21,Z23,Z31,Z32,Z14,Z24,Z34,Z41,Z42,Z43,c1,c2,c3,c4) << std::endl;
            std::cout << "Z12 = " << Z12 << '\n';


        }

        update(
            dt, &Z12, &Z12_n, &Z21, &Z21_n, &Z13, &Z13_n, &Z31, &Z31_n, &Z23, &Z23_n, 
            &Z32, &Z32_n, &Z41, &Z41_n, &Z42, &Z42_n, &Z43, &Z43_n, &Z14, &Z14_n, 
            &Z24, &Z24_n, &Z34, &Z34_n, &U12, &U13, &U21, &U23, &U31, &U32, &U14, &U24, &U34, &U41, &U42, &U43);
    
    }

    // Print the time it took to run this .cpp file.
    std::cout << "Finished" << std::time(nullptr)-start;
    return 0;
}