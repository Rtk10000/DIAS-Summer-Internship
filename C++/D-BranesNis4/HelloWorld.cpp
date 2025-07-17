#include <iostream>
#include <complex>
#include "eigen/Eigen/Dense"
typedef std::complex<double> R_or_C;

double acceleration(int j, double h)
{   
    double F;
    F += h;
    F += j;
    return F;
}
int main() {
    R_or_C john(12, 32);
    std::cout << john;
   
    return 0;
}


