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
const double sigma = .45;


// Repeat simulation for 1000 seconds.
const int simulation_repetitions = seconds_thermalised / delta_t;
// Number of D0-Branes
const int N = 4;
const int rows = N;
const int cols = N;

const long double g = 1/sqrt(N);

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
    R_or_C U = - 1.0/(4.0* g * g) * commutator_sum.trace().real();



    return std::abs(T + U);
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
            //double element = dist(rng);


            // For Complex elements Uncomment this
            long double real = dist(rng);
            long double imag = dist(rng);
            R_or_C element(real, imag);
            
            if (i == j)
            {
                Matrix(i,j) = real;
            }

            if (i != j) 
            {
                Matrix(i,j) = element;
                Matrix(j,i) = std::conj(element);
            }
        }
    }
    // Make matrix traceless by replacing last entry with the negative of the sum of the other entries.
    // Calculate the sum of diagonal elements in X1
    R_or_C diagonalSum = 0.0;
    for (int i = 0; i < rows-1; i++) 
    {
        diagonalSum += Matrix(i,i);
    }
    Matrix(rows - 1,cols - 1) = -diagonalSum;

    return Matrix;
}

int main()
{
    matrix zero_matrix = matrix::Zero(rows, cols);

    matrix V_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix,zero_matrix, zero_matrix, zero_matrix, zero_matrix};


    
    matrix X_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix,zero_matrix, zero_matrix, zero_matrix, zero_matrix};

    for (int i = 0; i < dim ; i++)
    {
        X_vector_new[i] = generateHermitianMatrix(rows, cols);

    }

    while (H(g, 
                    X_vector_new[0], X_vector_new[1], X_vector_new[2], X_vector_new[3], X_vector_new[4], X_vector_new[5], X_vector_new[6], X_vector_new[7], X_vector_new[8],
                    V_vector_new[0], V_vector_new[1], V_vector_new[2], V_vector_new[3], V_vector_new[4], V_vector_new[5], V_vector_new[6], V_vector_new[7], V_vector_new[8]) -1260> 0.01|| H(g, 
                    X_vector_new[0], X_vector_new[1], X_vector_new[2], X_vector_new[3], X_vector_new[4], X_vector_new[5], X_vector_new[6], X_vector_new[7], X_vector_new[8],
                    V_vector_new[0], V_vector_new[1], V_vector_new[2], V_vector_new[3], V_vector_new[4], V_vector_new[5], V_vector_new[6], V_vector_new[7], V_vector_new[8]) -1260 <-0.01)

    
    {
        for (int i = 0; i < dim ; i++)
        {
            X_vector_new[i] = generateHermitianMatrix(rows, cols);

        }

        std::cout << H(g, 
                            X_vector_new[0], X_vector_new[1], X_vector_new[2], X_vector_new[3], X_vector_new[4], X_vector_new[5], X_vector_new[6], X_vector_new[7], X_vector_new[8],
                          V_vector_new[0], V_vector_new[1], V_vector_new[2], V_vector_new[3], V_vector_new[4], V_vector_new[5], V_vector_new[6], V_vector_new[7], V_vector_new[8])<< '\n';

    }
    std::cout << "Done" << H(g, 
                            X_vector_new[0], X_vector_new[1], X_vector_new[2], X_vector_new[3], X_vector_new[4], X_vector_new[5], X_vector_new[6], X_vector_new[7], X_vector_new[8],
                          V_vector_new[0], V_vector_new[1], V_vector_new[2], V_vector_new[3], V_vector_new[4], V_vector_new[5], V_vector_new[6], V_vector_new[7], V_vector_new[8])<< '\n';

    
    for (matrix el : X_vector_new)
        {
            std::cout << el;

        }
    


    std:: fstream X_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D2-Branes/initial_X.txt", std:: ios:: out);
    X_vector_Export << std::fixed << std::setprecision(15);
    //Print to text file`
    for (matrix Matrix : X_vector_new)
        {
            X_vector_Export << Matrix << std::endl;
        }



 

    return 0;
}