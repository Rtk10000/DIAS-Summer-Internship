#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <random>
#include <complex>
#include <vector>
#include "eigen/Eigen/Dense"


int main()
{
    double F[9];
    F[0]= 9.0;
    F[1]= 9.0;
    F[2]= 9.0;
    F[3]= 9.0;
    F[4]= 9.0;
    F[5]= 9.0;
    F[6]= 2.0;
    F[7]= 9.0;
    F[8] = 9.0;

    double g[9];
    g[0]= 1.0;
    g[1]= 2.0;
    g[2]= 1.0;
    g[3]= 2.0;
    g[4]= 1.0;
    g[5]= 2.0;
    g[6]= 2.0;
    g[7]= 1.0;
    g[8]= 2.0;
    double e = F[8];

// Copy elements from X_vector_new to X_vector
    std::memcpy(F, g, sizeof(g)); 
    std::cout << e<< std::endl;

    for (double el : F)
    {
        std::cout << el;

    }

    return 0;
}