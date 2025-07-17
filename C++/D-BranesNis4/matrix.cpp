#include <iostream>
#include "eigen/Eigen/Dense"
const int N = 4;
typedef std::complex<double> R_or_C;
typedef Eigen:: Matrix<std::complex<double>, N, N> matrix;

int main()
{
    matrix t;
    t << 1,1,1,1,
         2,2,2,2,
         3,3,3,3,
         4,4,4,4;

    std::cout << t(3,1);




    return 0;
}