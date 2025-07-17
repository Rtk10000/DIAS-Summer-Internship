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
const double seconds_thermalised = 1000;
const double g = 0.01;

// Repeat simulation for 1000 seconds.
const int simulation_repetitions = seconds_thermalised / delta_t;
// Number of D0-Branes
const int N = 3;
const int rows = N;
const int cols = N;
const double sigma = 1;

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
    R_or_C T = 1.0/(2.0 ) * (V1 * V1 + V2 * V2 + V3 * V3 + V4 * V4 + V5 * V5 + V6 * V6 + V7 * V7 + V8 * V8 + V9 * V9).trace();

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
    R_or_C U = - g * g * 1.0/(4.0) * commutator_sum.trace();



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
    return g * g * commutator_sum;
}

// Acceleration of each coordinate matrix
matrix Acceleration2(const int i, matrix* X_vector, int rows, int cols, const double g)
{
    matrix commutator_sum = matrix::Zero(rows, cols);
    matrix temp_commutator = matrix::Zero(rows, cols);

    matrix X = X_vector[i];
    double c = X(0,1);
    double d = X(0,2);
    double e = X(1,2);
    double b = X(1,1);
    double a = X(0,0);
    for (int j = 0; j < dim; ++j)
    {
        if (i != j)
        {  
            matrix X_other = X_vector[j];
            double f = X_other(0,0);
            double g = X_other(1,1);
            double h = X_other(0,1);
            double i = X_other(0,2);
            double j = X_other(1,2);

            double zero_zero_entry = -2*a*h*h + 2*b*h*h + 4*e*i*h - 2*d*j*h + 2*c*f*h - 2*c*h*g - 4*a*i*i - 2*b*i*i - 2*c*i*j + 4*d*i*f + 2*d*i*g;
            double one_one_entry = 2*a*h*h - 2*b*h*h - 2*e*i*h + 4*d*j*h - 2*c*f*h + 2*c*h*g - 2*a*j*j - 4*b*j*j - 2*c*i*j + 2*e*j*f + 4*e*j*g;
            double zero_one_entry = -c*f*f + 3*d*j*f + a*f*h - b*f*h + 2*c*f*g - c*i*i - c*j*j - 3*a*i*j - 3*b*i*j + d*i*h + e*j*h - c*g*g + 3*e*i*g - a*h*g + b*h*g; 
            double zero_two_entry = -4*d*f*f + 4*a*i*f + 2*b*i*f + 3*c*j*f - 3*e*f*h - 4*d*f*g - d*j*j + e*i*j - d*h*h + c*i*h + 3*b*j*h - d*g*g - 3*e*h*g + 2*a*i*g + b*i*g; 
            double one_two_entry = -4*e*g*g + 3*c*i*g + 2*a*j*g + 4*b*j*g - 4*e*f*g - 3*d*h*g - e*i*i + d*i*j - e*f*f + a*j*f + 2*b*j*f - e*h*h + 3*a*i*h + c*j*h - 3*d*f*h;


            temp_commutator(0,0) = zero_zero_entry;
            temp_commutator(1,1) = one_one_entry;
            temp_commutator(0,1) = zero_one_entry;
            temp_commutator(0,2) = zero_two_entry;
            temp_commutator(1,2) = one_two_entry;
            temp_commutator(1,0) = temp_commutator(0,1);
            temp_commutator(2,0) = temp_commutator(0,2); 
            temp_commutator(2,1) = temp_commutator(1,2); 
 
            temp_commutator(2,2) = - zero_zero_entry - one_one_entry;  
            
            commutator_sum += temp_commutator;
        }   
    }
    return g * g * commutator_sum;
}

int main() 
{


    // Create  vectors to store the matrices
    matrix X_vector[dim];
    matrix V_vector[dim];
    matrix A_vector[dim];

    // Create a zero matrix in order to populate the V_vector with it.
    matrix zero_matrix = matrix::Zero(rows, cols);

/*
 //For testing reproducibility use these X values
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




    // Generate and store V1, V2, V3, V4, V5, V6, V7, V8, and V9
    for (int i = 0; i < dim; ++i) 
    {
        V_vector[i] = zero_matrix;
    }  

    // Generate and store A1, A2, A3, A4, A5, A6, A7, A8, and A9
    for (int i = 0; i < dim; ++i) 
    {
        A_vector[i] = Acceleration2(i, X_vector, rows, cols, g);
    }  


    // Create  vectors to store the new matrices
    matrix X_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix V_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix A_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};

    // Write simulation to thermalise system
    for (int j = 0; j < seconds_thermalised / delta_t; ++j)
    {

        // velocity Verlet 1 to get new positions from old positions, momentums and rate of change of momentums

        for (int i = 0; i < 9; ++i)
        {
            X_vector_new[i] = X_vector[i] + V_vector[i] * delta_t + 0.5 * A_vector[i] * delta_t * delta_t;
        }

        // Generate and store new A1, A2, A3, A4, A5, A6, A7, A8, and A9
        for (int i = 0; i < 9; ++i) 
        {
            A_vector_new[i] = Acceleration2( i, X_vector_new, rows, cols, g);
        }  
        
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


        }

    }




       // Export initial X/V/A_vector to text files to be analysed in python.
    std:: fstream X2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D0-BranesNis3/thermalised_XNis3.txt", std:: ios:: out);
    X2_vector_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (matrix Matrix : X_vector_new)
    {
        X2_vector_Export << Matrix << std::endl;
    }


    std:: fstream V2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D0-BranesNis3/thermalised_VNis3.txt", std:: ios:: out);
    V2_vector_Export << std::fixed << std::setprecision(15);
    for (matrix Matrix : V_vector_new)
    {
        V2_vector_Export << Matrix << std::endl;
    }

    std:: fstream A2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D0-BranesNis3/thermalised_ANis3.txt", std:: ios:: out);
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