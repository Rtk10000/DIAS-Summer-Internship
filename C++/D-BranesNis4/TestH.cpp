#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <random>
#include <complex>
#include <vector>
#include "eigen/Eigen/Dense"
double start = std::time(nullptr);

// Define timestep
const long double delta_t = 1e-4;
const long double seconds_thermalised = 1000;

// Repeat simulation for 1000 seconds.
const int simulation_repetitions = seconds_thermalised / delta_t;
// Number of D0-Branes
const int N = 4;
const long double g = 1.0/sqrt(N);
const int rows = N;
const int cols = N;

// Dimension of space
const int dim = 9;


//typedef double R_or_C;
//typedef Eigen:: Matrix<double, N, N> matrix;


//To go from real matrices to complex matrices delete the above typedefs and use these typedefs.
//And change the commented out code in generateHermitianMatrix()

typedef std::complex<long double> R_or_C;
typedef Eigen:: Matrix<std::complex<long double>, N, N> matrix;


matrix commutator(matrix A, matrix B)
{
    return A * B - B * A;
} 


// Cillian's Hamiltonian
long double H(
    long double g, 
    matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7, matrix X8, matrix X9,
    matrix V1, matrix V2, matrix V3, matrix V4, matrix V5, matrix V6, matrix V7, matrix V8, matrix V9)
{
    // Compute kinetic energy T
    R_or_C T = 1.0/(2.0 * g * g) * (V1 * V1 + V2 * V2 + V3 * V3 + V4 * V4 + V5 * V5 + V6 * V6 + V7 * V7 + V8 * V8 + V9 * V9).trace();

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
    R_or_C U = - 1.0/(4.0 * g * g) * commutator_sum.trace();


    return std::abs(T + U.real());
}


int main() 
{
    std::cout << "g" << g << '\n';
    // Create  vectors to store the matrices
    matrix X_vector[dim] = {matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols)};
    matrix V_vector[dim]= {matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols)};
    matrix A_vector[dim]= {matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols)};

    // Create a zero matrix in order to populate the V_vector with it.
    matrix zero_matrix = matrix::Zero(rows, cols);


// For testing reproducibility use these X values
    std::ifstream inputX("TestX.txt");
    if (!inputX.is_open()) {
        std::cerr << "Failed to open X initial file." << std::endl;
        return 1;
    }

    // Read the values from the file and store them in the matrices
    for (int i = 0; i < dim; ++i) 
    {

        for (int row = 0; row < rows; ++row) 
        {
            for (int col = 0; col < cols; ++col) 
            {
                 inputX >> X_vector[i](row, col);
            }
        }
    }

    // Close the input file
    inputX.close();


    // For testing reproducibility use these X values
    std::ifstream inputV("TestV.txt");
    if (!inputV.is_open()) {
        std::cerr << "Failed to open X initial file." << std::endl;
        return 1;
    }

    // Read the values from the file and store them in the matrices
    for (int i = 0; i < dim; ++i) 
    {

        for (int row = 0; row < rows; ++row) 
        {
            for (int col = 0; col < cols; ++col) 
            {
                 inputV >> V_vector[i](row, col);
            }
        }
    }

    // Close the input file
    inputV.close();


    for(matrix el : X_vector)
    {
        std::cout << '\n' <<el;
    }



    std::cout << "H" << std::setprecision(15) << H(g, 
                            X_vector[0], X_vector[1], X_vector[2], X_vector[3], X_vector[4], X_vector[5], X_vector[6], X_vector[7], X_vector[8],
                          V_vector[0], V_vector[1], V_vector[2], V_vector[3], V_vector[4], V_vector[5], V_vector[6], V_vector[7], V_vector[8]);


    return 0;
}