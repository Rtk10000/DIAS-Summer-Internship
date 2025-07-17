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
const double delta_t = 1e-5;
const double seconds_thermalised = 10;
const double g = 1;

// Repeat simulation for 1000 seconds.
const int simulation_repetitions = seconds_thermalised / delta_t;
// Number of D0-Branes
const int N = 2;
const int rows = N;
const int cols = N;

// Dimension of space
const int dim = 9;


typedef double R_or_C;
typedef Eigen::Matrix<double, N, N> matrix;


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
    R_or_C T = 1.0 / (2.0 * pow(g, 2)) * (V1 * V1 + V2 * V2 + V3 * V3 + V4 * V4 + V5 * V5 + V6 * V6 + V7 * V7 + V8 * V8 + V9 * V9).trace();

    matrix X[9] = { X1,X2,X3,X4,X5,X6,X7,X8,X9 };

    matrix commutator_sum = matrix::Zero(rows, cols);
    for (int i = 0; i < 9; i++)
    {
        for (int j = 0; j < 9; j++)
        {
            if (i != j)
            {
                commutator_sum += commutator(X[i], X[j]) * commutator(X[i], X[j]); //can likely be more efficient by less function calls

            }
        }
    }
    R_or_C U = -1.0 / (4.0 * pow(g, 2)) * commutator_sum.trace();



    return std::abs(T + U);
}

// Cillian's Gauss' law
matrix gauss_law(
    matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7, matrix X8, matrix X9,
    matrix V1, matrix V2, matrix V3, matrix V4, matrix V5, matrix V6, matrix V7, matrix V8, matrix V9)
{
    matrix result;
    for (int i = 0; i < 9; i++)
    {
        result = commutator(X1, V1) + commutator(X2, V2) + commutator(X3, V3) + commutator(X4, V4) + commutator(X5, V5) +
            commutator(X6, V6) + commutator(X7, V7) + commutator(X8, V8) + commutator(X9, V9);
    }
    return result;
}


// Random Hermitian matrix

static std::random_device rd;

static std::mt19937 rng(std::time(nullptr));
std::normal_distribution<double> dist(0.0, 1);

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

            Matrix(i, j) = element;
            if (i != j)
            {
                Matrix(j, i) = element;
            }
        }
    }
    // Make matrix traceless by replacing last entry with the negative of the sum of the other entries.
    // Calculate the sum of diagonal elements in X1
    R_or_C diagonalSum = 0.0;
    for (int i = 0; i < rows - 1; ++i)
    {
        diagonalSum += Matrix(i, i);
    }
    Matrix(rows - 1, cols - 1) = -diagonalSum;

    return Matrix;
}


// Acceleration of each coordinate matrix
matrix Acceleration(const int j, matrix* X_vector, int rows, int cols, const double g)
{
    matrix commutator_sum = matrix::Zero(rows, cols);

    matrix X = X_vector[j];
    for (int i = 0; i < dim; ++i)
    {
        if (i != j)
        {
            matrix temp_commutator = commutator(X_vector[i], commutator(X, X_vector[i]));


            commutator_sum += temp_commutator;

        }
        // Comment out the line below to go back to the lagrangian e.q.m. With no -1/(g^2)
        // commutator_sum = -1.0 / (g * g) * commutator_sum;
    }
    return commutator_sum;
}

// Acceleration of each coordinate matrix
matrix Acceleration2(const double aX, const double bX, const double a1, const double b1, const double a2, const double b2, 
    const double a3, const double b3, const double a4, const double b4, const double a5, const double b5, const double a6, const double b6, 
    const double a7, const double b7, const double a8, const double b8, const double g)
{
    matrix commutator_sum = matrix::Zero(rows, cols);
    matrix temp_commutator = matrix::Zero(rows, cols);

    double matrix_entries[16] = { a1, b1, a2,   b2, a3,   b3,   a4,   b4,   a5,   b5,   a6,  b6, a7,  b7,  a8,  b8 };

    for (int j = 0; j < dim-1; ++j)
    {

        temp_commutator(0, 0) = -4 * matrix_entries[j + 1] * matrix_entries[j + 1] * aX + 4 * matrix_entries[j + 1] * matrix_entries[j] * bX;
        temp_commutator(0, 1) = 4 * matrix_entries[j] * (matrix_entries[j + 1] * aX - matrix_entries[j] * bX);
        temp_commutator(1, 0) = temp_commutator(0, 1);
        temp_commutator(1, 1) = -temp_commutator(0, 0);

        commutator_sum += temp_commutator;
      
    }
    return commutator_sum;
}

int main()
{


    // Create  vectors to store the matrices
    matrix X_vector[dim];
    matrix V_vector[dim];
    matrix A_vector[dim];

    /*
    // For testing reproducibility use these X values
        std::ifstream inputX("initial_X.txt");
        if (!inputX.is_open()) {
            std::cerr << "Failed to open the file." << std::endl;
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
    */

    // Generate and store X1, X2, X3, X4, X5, X6, X7, X8, and X9
    for (int i = 0; i < dim; ++i)
    {
        X_vector[i] = generateHermitianMatrix(rows, cols);
    }



    // Create a zero matrix in order to populate the V_vector with it.
    matrix zero_matrix = matrix::Zero(rows, cols);

    // Generate and store V1, V2, V3, V4, V5, V6, V7, V8, and V9
    for (int i = 0; i < dim; ++i)
    {
        V_vector[i] = zero_matrix;
    }

    // Generate and store A1, A2, A3, A4, A5, A6, A7, A8, and A9
    double a0 = X_vector[0](0, 0);
    double b0 = X_vector[0](0, 1);
    double a1 = X_vector[1](0, 0);
    double b1 = X_vector[1](0, 1);
    double a2 = X_vector[2](0, 0);
    double b2 = X_vector[2](0, 1);
    double a3 = X_vector[3](0, 0);
    double b3 = X_vector[3](0, 1);
    double a4 = X_vector[4](0, 0);
    double b4 = X_vector[4](0, 1);
    double a5 = X_vector[5](0, 0);
    double b5 = X_vector[5](0, 1);
    double a6 = X_vector[6](0, 0);
    double b6 = X_vector[6](0, 1);
    double a7 = X_vector[7](0, 0);
    double b7 = X_vector[7](0, 1);
    double a8 = X_vector[8](0, 0);
    double b8 = X_vector[8](0, 1);
    // Generate and store new A1, A2, A3, A4, A5, A6, A7, A8, and A9

    A_vector[0] = Acceleration2(a0, b0, a1, b1, a2, b2, a3, b3, a4, b4, a5, b5, a6, b6, a7, b7, a8, b8, g);
    A_vector[1] = Acceleration2(a1, b1, a0, b0, a2, b2, a3, b3, a4, b4, a5, b5, a6, b6, a7, b7, a8, b8, g);
    A_vector[2] = Acceleration2(a2, b2, a1, b1, a0, b0, a3, b3, a4, b4, a5, b5, a6, b6, a7, b7, a8, b8, g);
    A_vector[3] = Acceleration2(a3, b3, a1, b1, a2, b2, a0, b0, a4, b4, a5, b5, a6, b6, a7, b7, a8, b8, g);
    A_vector[4] = Acceleration2(a4, b4, a1, b1, a2, b2, a3, b3, a0, b0, a5, b5, a6, b6, a7, b7, a8, b8, g);
    A_vector[5] = Acceleration2(a5, b5, a1, b1, a2, b2, a3, b3, a4, b4, a0, b0, a6, b6, a7, b7, a8, b8, g);
    A_vector[6] = Acceleration2(a6, b6, a1, b1, a2, b2, a3, b3, a4, b4, a5, b5, a0, b0, a7, b7, a8, b8, g);
    A_vector[7] = Acceleration2(a7, b7, a1, b1, a2, b2, a3, b3, a4, b4, a5, b5, a6, b6, a0, b0, a8, b8, g);
    A_vector[8] = Acceleration2(a8, b8, a1, b1, a2, b2, a3, b3, a4, b4, a5, b5, a6, b6, a7, b7, a0, b0, g);



    // Export initial X/V/A_vector to text files to be analysed in python.

    std::fstream X_vector_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/X_vector_test4.txt", std::ios::out);
    X_vector_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (matrix Matrix : X_vector)
    {
        X_vector_Export << Matrix << std::endl;
    }

    std::fstream V_vector_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/V_vector_test4.txt", std::ios::out);
    V_vector_Export << std::fixed << std::setprecision(15);
    for (matrix Matrix : V_vector)
    {
        V_vector_Export << Matrix << std::endl;
    }

    std::fstream A_vector_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/A_vector_test4.txt", std::ios::out);
    A_vector_Export << std::fixed << std::setprecision(15);
    // Print to text file
    for (matrix Matrix : A_vector)
    {
        A_vector_Export << Matrix << std::endl;
    }

    // Create  vectors to store the new matrices
    matrix X_vector_new[9] = { zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix };
    matrix V_vector_new[9] = { zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix };
    matrix A_vector_new[9] = { zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix };

    // Write simulation to thermalise system
    for (int j = 0; j < seconds_thermalised / delta_t; ++j)
    {

        // velocity Verlet 1 to get new positions from old positions, momentums and rate of change of momentums

        for (int i = 0; i < 9; ++i)
        {
            X_vector_new[i] = X_vector[i] + V_vector[i] * delta_t + 0.5 * A_vector[i] * delta_t * delta_t;
        }

        double a0 = X_vector_new[0](0, 0);
        double b0 = X_vector_new[0](0, 1);
        double a1 = X_vector_new[1](0, 0);
        double b1 = X_vector_new[1](0, 1);
        double a2 = X_vector_new[2](0, 0);
        double b2 = X_vector_new[2](0, 1);
        double a3 = X_vector_new[3](0, 0);
        double b3 = X_vector_new[3](0, 1);
        double a4 = X_vector_new[4](0, 0);
        double b4 = X_vector_new[4](0, 1);
        double a5 = X_vector_new[5](0, 0);
        double b5 = X_vector_new[5](0, 1);
        double a6 = X_vector_new[6](0, 0);
        double b6 = X_vector_new[6](0, 1);
        double a7 = X_vector_new[7](0, 0);
        double b7 = X_vector_new[7](0, 1);
        double a8 = X_vector_new[8](0, 0);
        double b8 = X_vector_new[8](0, 1);
        // Generate and store new A1, A2, A3, A4, A5, A6, A7, A8, and A9

        A_vector_new[0] = Acceleration2( a0, b0, a1, b1, a2, b2, a3,  b3,  a4, b4, a5,  b5,  a6,  b6, a7,  b7,  a8,  b8,  g);
        A_vector_new[1] = Acceleration2(a1, b1, a0, b0, a2, b2, a3, b3, a4, b4, a5, b5, a6, b6, a7, b7, a8, b8, g);
        A_vector_new[2] = Acceleration2(a2, b2, a1, b1, a0, b0, a3, b3, a4, b4, a5, b5, a6, b6, a7, b7, a8, b8, g);
        A_vector_new[3] = Acceleration2(a3, b3, a1, b1, a2, b2, a0, b0, a4, b4, a5, b5, a6, b6, a7, b7, a8, b8, g);
        A_vector_new[4] = Acceleration2(a4, b4, a1, b1, a2, b2, a3, b3, a0, b0, a5, b5, a6, b6, a7, b7, a8, b8, g);
        A_vector_new[5] = Acceleration2(a5, b5, a1, b1, a2, b2, a3, b3, a4, b4, a0, b0, a6, b6, a7, b7, a8, b8, g);
        A_vector_new[6] = Acceleration2(a6, b6, a1, b1, a2, b2, a3, b3, a4, b4, a5, b5, a0, b0, a7, b7, a8, b8, g);
        A_vector_new[7] = Acceleration2(a7, b7, a1, b1, a2, b2,  a3, b3, a4, b4, a5, b5, a6, b6, a0, b0, a8, b8, g);
        A_vector_new[8] = Acceleration2(a8, b8, a1, b1, a2, b2, a3, b3, a4, b4, a5, b5, a6, b6, a7, b7, a0, b0, g);
   

        // Use Velocity Verlet 2 to get new momentums
        for (int i = 0; i < 9; ++i)
        {
            V_vector_new[i] = V_vector[i] + 0.5 * (A_vector_new[i] + A_vector[i]) * delta_t;
        }

        // Copy elements from X_vector_new to X_vector
        std::memcpy(X_vector, X_vector_new, sizeof(X_vector_new));

        // Copy elements from X_vector_new to X_vector
        std::memcpy(V_vector, V_vector_new, sizeof(V_vector_new));

        // Copy elements from X_vector_new to X_vector
        std::memcpy(A_vector, A_vector_new, sizeof(A_vector_new));

        if (j % 10000 == 0)
        {
            //for (matrix el : V_vector)
            //{
            //    std::cout <<"Ideal " << el << std::endl;
            //}
            //std::cout  << std::endl << gauss_law(X_vector_new[0], X_vector_new[1], X_vector_new[2], X_vector_new[3], X_vector_new[4], X_vector_new[5], X_vector_new[6], X_vector_new[7], X_vector_new[8],
            //                       V_vector_new[0], V_vector_new[1], V_vector_new[2], V_vector_new[3], V_vector_new[4], V_vector_new[5], V_vector_new[6], V_vector_new[7], V_vector_new[8]);

            std::cout << std::endl;
            std::cout << "H" << std::setprecision(15) << H(1.0,
                X_vector_new[0], X_vector_new[1], X_vector_new[2], X_vector_new[3], X_vector_new[4], X_vector_new[5], X_vector_new[6], X_vector_new[7], X_vector_new[8],
                V_vector_new[0], V_vector_new[1], V_vector_new[2], V_vector_new[3], V_vector_new[4], V_vector_new[5], V_vector_new[6], V_vector_new[7], V_vector_new[8]);

            /*
                        for ( matrix Matrix : X_vector_new)
                        {
                            X_vector_Export << Matrix << std::endl;
                        }

                        for ( matrix Matrix : V_vector_new)
                        {
                            V_vector_Export << Matrix << std::endl;
                        }

                        for ( matrix Matrix : A_vector_new)
                        {
                            A_vector_Export << Matrix << std::endl;
                        }
            */
        }

    }


    X_vector_Export.close();
    V_vector_Export.close();
    A_vector_Export.close();

    // Export initial X/V/A_vector to text files to be analysed in python.
    std::fstream X2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D0-Branes/thermalised_X.txt", std::ios::out);
    X2_vector_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (matrix Matrix : X_vector_new)
    {
        X2_vector_Export << Matrix << std::endl;
    }


    std::fstream V2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D0-Branes/thermalised_V.txt", std::ios::out);
    V2_vector_Export << std::fixed << std::setprecision(15);
    for (matrix Matrix : V_vector_new)
    {
        V2_vector_Export << Matrix << std::endl;
    }

    std::fstream A2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D0-Branes/thermalised_A.txt", std::ios::out);
    A2_vector_Export << std::fixed << std::setprecision(15);
    // Print to text file
    for (matrix Matrix : A_vector_new)
    {
        A2_vector_Export << Matrix << std::endl;
    }


    X2_vector_Export.close();
    V2_vector_Export.close();
    A2_vector_Export.close();



    return 0;
}