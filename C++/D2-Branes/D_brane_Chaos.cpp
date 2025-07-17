#include <iostream>
#include <random>
#include <ctime>
#include <string>
#include <complex>
#include <fstream>
#include <sstream>
#include "eigen/Eigen/Dense"

const int N = 3;
typedef Eigen:: Matrix<std:: complex<double>, N, N> matrix;
typedef std:: complex<double> complex;

/*
Reserving memory for my P's, Q's, X's etc. So the P's will be a row/column vector where all components are real,
and the same will be true for the Q's. Then the Z's can be constructed from P + complex(0,1)*Q. there should be 6
scalar Z's and 6 vector Z's. Then I will have 27 scalar X values which should be easy to implement just as complex(a,b) for each
and then I will have the 9 NxN X matrices. These can just be random complex matrices to start. All degrees of freedom can be random
to start and the dotted coordinates can be zero like the case of D0 branes. Stick with this until Swapno confirms this satisfies the new 
Gauss Law constraint. Also must ask swapno how the real and complex derivatives can be combined to give rise to the full equation of motion
of the complex valued degrees of freedom. 
N.B. The X matrices should be hermitian but they need not be traceless in this case.  
*/

std:: mt19937 rng(std::time(nullptr));
std:: normal_distribution<long double> gauss_dist(0, 1);

// Initialise everything as a random variable for now until swapno gets back and says what things should be
// Scalar Z's 
double P12 = gauss_dist(rng), P13 = gauss_dist(rng), P21 = gauss_dist(rng), P23 = gauss_dist(rng), P31 = gauss_dist(rng), P32 = gauss_dist(rng);
double Q12 = gauss_dist(rng), Q13 = gauss_dist(rng), Q21 = gauss_dist(rng), Q23 = gauss_dist(rng), Q31 = gauss_dist(rng), Q32 = gauss_dist(rng);
complex Z12 = complex(P12,Q12);
complex Z13 = complex(P13,Q13);
complex Z21 = complex(P21,Q21);
complex Z23 = complex(P23,Q23); 
complex Z31 = complex(P31,Q31); 
complex Z32 = complex(P32,Q32); 

// Vector Z's
Eigen::Matrix<double, N, 1> P41 = Eigen::Matrix<double, N, 1>::Random();
Eigen::Matrix<double, N, 1> P42 = Eigen::Matrix<double, N, 1>::Random();
Eigen::Matrix<double, N, 1> P43 = Eigen::Matrix<double, N, 1>::Random();
Eigen::Matrix<double, N, 1> Q41 = Eigen::Matrix<double, N, 1>::Random();
Eigen::Matrix<double, N, 1> Q42 = Eigen::Matrix<double, N, 1>::Random(); 
Eigen::Matrix<double, N, 1> Q43 = Eigen::Matrix<double, N, 1>::Random();

Eigen::Matrix<double, 1, N> P14 = Eigen::Matrix<double, 1, N>::Random();
Eigen::Matrix<double, 1, N> P24 = Eigen::Matrix<double, 1, N>::Random(); 
Eigen::Matrix<double, 1, N> P34 = Eigen::Matrix<double, 1, N>::Random();
Eigen::Matrix<double, 1, N> Q14 = Eigen::Matrix<double, 1, N>::Random();
Eigen::Matrix<double, 1, N> Q24 = Eigen::Matrix<double, 1, N>::Random();
Eigen::Matrix<double, 1, N> Q34 = Eigen::Matrix<double, 1, N>::Random();

Eigen::Matrix<std::complex<double>, N, 1> Z41 = P41 + complex(0,1)*Q41;
Eigen::Matrix<std::complex<double>, N, 1> Z42 = P42 + complex(0,1)*Q42;
Eigen::Matrix<std::complex<double>, N, 1> Z43 = P43 + complex(0,1)*Q43;
Eigen::Matrix<std::complex<double>, 1, N> Z14 = P14 + complex(0,1)*Q14;
Eigen::Matrix<std::complex<double>, 1, N> Z24 = P24 + complex(0,1)*Q24;
Eigen::Matrix<std::complex<double>, 1, N> Z34 = P34 + complex(0,1)*Q34;

// Scalar X's
complex X11, X12, X13, X14, X15, X16, X17, X18, X19;
complex X21, X22, X23, X24, X25, X26, X27, X28, X29;
complex X31, X32, X33, X34, X35, X36, X37, X38, X39;

// Matrix X's
matrix X41, X42, X43, X44, X45, X46, X47, X48, X49;

int main()
{
    return 0;
}