#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <random>
#include <complex>
#include <vector>
#include "eigen/Eigen/Dense"

const int N = 4;

typedef std::complex<long double> R_or_C;
typedef std::complex<long double> Complex;

typedef Eigen:: Matrix<std::complex<long double>, N, 1> col_vector;
typedef Eigen:: Matrix<std::complex<long double>, 1, N> row_vector;
typedef Eigen:: Matrix<std::complex<long double>, N, N> matrix;

long double delta_t = 0.1;
// Initialise a random engine to be used when creating random initial Z coordinates 
std:: mt19937 rng(std::time(nullptr));
std:: normal_distribution<long double> gauss_dist(0, 1);

R_or_C I(0,1);
int main()
{

     col_vector Z41;
    Z41(0,0) = Complex(1,1);
    Z41(1,0) = Complex(1,1);
    Z41(2,0) = Complex(1,1);
    Z41(3,0) = Complex(1,1);

    row_vector Z34;
    Z34(0,0) = Complex(1,3);
    Z34(0,1) = Complex(2,3);
    Z34(0,2) = Complex(3,3);
    Z34(0,3) = Complex(4,3);


    std::cout << Z41 << '\n';




    return ;
}