// This code is for creating reproducible random initial conditions into initial_X.txt


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
const double seconds_thermalised = 20;
const double g = 1;
const double sigma = 0.1;

// Repeat simulation for 1000 seconds.
const int simulation_repetitions = seconds_thermalised / delta_t;
// Number of D0-Branes
const int N = 2;
const int rows = N;
const int cols = N;

// Dimension of space
const int dim = 9;


typedef double R_or_C;
typedef Eigen:: Matrix<double, N, N> matrix;


//To go from real matrices to complex matrices delete the above typedefs and use these typedefs.
//And change the commented out code in generateHermitianMatrix()

//typedef std::complex<double> R_or_C;
//typedef Eigen:: Matrix<std::complex<double>, N, N> matrix;




matrix commutator(matrix A, matrix B)
{
    return A * B - B * A;
} 


// Cillian's Hamiltonian
double H(
    double g, 
    matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7, matrix X8, matrix X9,
    matrix V1, matrix V2, matrix V3, matrix V4, matrix V5, matrix V6, matrix V7, matrix V8, matrix V9)
{
    // Compute kinetic energy T
    R_or_C T = 1.0/(2.0 * pow(g,2)) * (V1 * V1 + V2 * V2 + V3 * V3 + V4 * V4 + V5 * V5 + V6 * V6 + V7 * V7 + V8 * V8 + V9 * V9).trace();

    matrix X[9] = {X1,X2,X3,X4,X5,X6,X7,X8,X9}; 

    matrix commutator_sum = matrix::Zero(rows, cols);  
    for (int i = 0; i < 9; i++)
    {
        for (int j = 0; j < 9; j++)
        {
            if(i != j)
            {
                commutator_sum += commutator(X[i],X[j])*commutator(X[i],X[j]); //can likely be more efficient by less function calls

            }
        }
    }
    R_or_C U = - 1.0/(4.0*pow(g,2)) * commutator_sum.trace();



    return std:: abs(T + U);
}

// Cillian's Gauss' law
matrix gauss_law(
    matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7, matrix X8, matrix X9,
    matrix V1, matrix V2, matrix V3, matrix V4, matrix V5, matrix V6, matrix V7, matrix V8, matrix V9)
{
    matrix result;
    for (int i = 0; i < 9; i ++)
    {
        result = commutator(X1,V1) + commutator(X2,V2) + commutator(X3,V3) + commutator(X4,V4) + commutator(X5,V5) +
        commutator(X6,V6) + commutator(X7,V7) + commutator(X8,V8) + commutator(X9,V9);
    }
    return result;
}


// Random Hermitian matrix

static std::random_device rd;

static std::mt19937 rng(std::time(nullptr)); 
std::normal_distribution<double> dist(0.0, sigma);

matrix generateHermitianMatrix(int rows, int cols) 
{   
    matrix Matrix(rows, cols);

    for (int i = 0; i < rows; ++i) 
    {
        for (int j = 0; j <= i; ++j) 
        {
            double element = dist(rng);
            //double real = dist(rng);
            //double imag = dist(rng);
            //R_or_C element(real, imag);

            /* For Complex elements Uncomment this
            double real = dist(rng);
            double imag = dist(rng);
            R_or_C element(real, imag);
            */

            Matrix(i,j) = element;
            if (i != j) 
            {
                Matrix(j,i) = std::conj(element);
            }
        }
    }
    // Make matrix traceless by replacing last entry with the negative of the sum of the other entries.
    // Calculate the sum of diagonal elements in X1
    R_or_C diagonalSum = 0.0;
    for (int i = 0; i < rows-1; ++i) 
    {
        diagonalSum += Matrix(i,i);
    }
    Matrix(rows - 1,cols - 1) = -diagonalSum;

    return Matrix;
}

int main()
{
    matrix X_vector[dim];
    // Generate and store X1, X2, X3, X4, X5, X6, X7, X8, and X9
    for (int i = 0; i < dim; ++i) 
    {
        X_vector[i] = generateHermitianMatrix(rows, cols);
    }

    std:: fstream X_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D0-BranesNis2/initial_X.txt", std:: ios:: out);
    X_vector_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (matrix Matrix : X_vector)
        {
            X_vector_Export << Matrix << std::endl;
        }


    return 0;
}