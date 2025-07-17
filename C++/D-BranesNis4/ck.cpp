#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <random>
#include <complex>
#include <vector>
#include "eigen/Eigen/Dense"

// Define timestep
const double delta_t = 1e-4;
const double seconds_perturbed = 1;

// Repeat simulation for 1000 seconds.
const int simulation_repetitions = seconds_perturbed / delta_t;
// Number of D0-Branes
const int N = 4;
const int rows = N;
const int cols = N;

const double g = 1/sqrt(N);

// Dimension of space
const int dim = 9;

typedef std::complex<double> R_or_C;
typedef Eigen:: Matrix<std::complex<double>, N, N> matrix;

//typedef double R_or_C;
//typedef Eigen::Matrix2d matrix;

matrix commutator(matrix A, matrix B)
{
    return A * B - B * A;
} 

matrix anti_commutator(matrix A, matrix B)
{
    return A * B + B * A;
} 

static std::random_device rd;
static std::mt19937 rng(std::time(nullptr)); 
std::normal_distribution<double> dist(0.0, 1e-8);

int main()
{

    double c_1 = dist(rng);
    double c_2 = dist(rng);
    std::cout << c_1 << " " << c_2;

    return 0;
}