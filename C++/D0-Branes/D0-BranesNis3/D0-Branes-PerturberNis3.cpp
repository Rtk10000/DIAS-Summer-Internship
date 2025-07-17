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
const double g = 1;

// Repeat simulation for 1000 seconds.
const int simulation_repetitions = seconds_perturbed / delta_t;
// Number of D0-Branes
const int N = 3;
const int rows = N;
const int cols = N;

// Dimension of space
const int dim = 9;

//typedef std::complex<double> R_or_C;
//typedef Eigen:: Matrix<std::complex<double>, N, N> matrix;

typedef double R_or_C;
typedef Eigen:: Matrix<double, N, N> matrix;

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

// Cillian's Hamiltonian
double H(
    double g, 
    matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7, matrix X8, matrix X9,
    matrix V1, matrix V2, matrix V3, matrix V4, matrix V5, matrix V6, matrix V7, matrix V8, matrix V9)
{
    // Compute kinetic energy T
    R_or_C T = 1.0/(2.0) * (V1*V1 + V2*V2 + V3*V3 + V4*V4 + V5*V5 + V6*V6 + V7*V7 + V8*V8 + V9*V9).trace();

    matrix X[9] = {X1,X2,X3,X4,X5,X6,X7,X8,X9}; 

    matrix commutator_sum = matrix::Zero(rows, cols);  
    for (int i = 0; i < 9; i++)
    {
        for (int j = 0; j < 9; j++)
        {
            if(i == j)
                continue;
            commutator_sum += commutator(X[i],X[j])*commutator(X[i],X[j]); //can likely be more efficient by less function calls
        }
    }
    R_or_C U = - g * g * 1.0/(4.0) * commutator_sum.trace();
    return std:: abs(T + U);
}

// Cillian's Gauss' law
matrix gauss_law(matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7, matrix X8, matrix X9,
                 matrix V1, matrix V2, matrix V3, matrix V4, matrix V5, matrix V6, matrix V7, matrix V8, matrix V9)
{
    matrix result;
    for (int i = 0; i < 9; i++)
    {
        result = commutator(X1,V1) + commutator(X2,V2) + commutator(X3,V3) + commutator(X4,V4) + commutator(X5,V5) +
        commutator(X6,V6) + commutator(X7,V7) + commutator(X8,V8) + commutator(X9,V9);
    }
    return result;
}

// Acceleration of each coordinate matrix
matrix PerturbingAcceleration(const int i, matrix* X_vector, int rows, int cols, const double g, const double c_1, const double c_2)
{
    matrix commutator_sum = matrix::Zero(rows, cols);
    matrix sum_X = matrix::Zero(rows, cols);
    matrix temp_commutator = matrix::Zero(rows, cols);

    matrix X = X_vector[i];
    double c = X(0,1);
    double d = X(0,2);
    double e = X(1,2);
    double b = X(1,1);
    double a = X(0,0);

    // The below 6 lines are doing the perturbing.
    for( int k = 0; k < dim; k++)
    {
        sum_X += X_vector[k] * X_vector[k];
    }
    matrix perturbation = 2.0 * c_1 * X + 2.0 * g * g *  c_2 * anti_commutator(X, sum_X );
    commutator_sum += perturbation;

    for (int j = 0; j < dim; j++)
    {

        sum_X = matrix::Zero(rows, cols);

        if (j != i)
        {  
            matrix X_other = X_vector[j];
            double f = X_other(0,0);
            double g = X_other(1,1);
            double h = X_other(0,1);
            double i = X_other(0,2);
            double j = X_other(1,2);

            temp_commutator(0,0) = -2*a*h*h + 2*b*h*h + 4*e*i*h - 2*d*j*h + 2*c*f*h - 2*c*h*g - 4*a*i*i - 2*b*i*i - 2*c*i*j + 4*d*i*f + 2*d*i*g;
            temp_commutator(1,1) = 2*a*h*h - 2*b*h*h - 2*e*i*h + 4*d*j*h - 2*c*f*h + 2*c*h*g - 2*a*j*j - 4*b*j*j - 2*c*i*j + 2*e*j*f + 4*e*j*g;
            temp_commutator(0,1) = -c*f*f + 3*d*j*f + a*f*h - b*f*h + 2*c*f*g - c*i*i - c*j*j - 3*a*i*j - 3*b*i*j + d*i*h + e*j*h - c*g*g + 3*e*i*g - a*h*g + b*h*g;
            temp_commutator(0,2) = -4*d*f*f + 4*a*i*f + 2*b*i*f + 3*c*j*f - 3*e*f*h - 4*d*f*g - d*j*j + e*i*j - d*h*h + c*i*h + 3*b*j*h - d*g*g - 3*e*h*g + 2*a*i*g + b*i*g;
            temp_commutator(1,2) = -4*e*g*g + 3*c*i*g + 2*a*j*g + 4*b*j*g - 4*e*f*g - 3*d*h*g - e*i*i + d*i*j - e*f*f + a*j*f + 2*b*j*f - e*h*h + 3*a*i*h + c*j*h - 3*d*f*h;
            temp_commutator(1,0) = temp_commutator(0,1);
            temp_commutator(2,0) = temp_commutator(0,2); 
            temp_commutator(2,1) = temp_commutator(1,2); 
 
            temp_commutator(2,2) = - temp_commutator(0,0) - temp_commutator(1,1);  

            commutator_sum += g * g * temp_commutator;
        }   

    // Comment out the line below to go back to the lagrangian e.q.m. With no -1/(g^2)
    //commutator_sum = -1.0 / (g * g) * commutator_sum;
    }
    return commutator_sum;
}


int main() 
{


    R_or_C c_1 = dist(rng);
    R_or_C c_2 = dist(rng);
    // Create  vectors to store the matrices
    matrix X2_vector[dim];
    matrix V2_vector[dim];
    matrix A2_vector[dim];



    // Generate and store X1, X2, X3, X4, X5, X6, X7, X8, and X9
        // Create an array to store the matrices
    // Open the text file for reading
    std::ifstream inputX("Thermalised_XNis3.txt");
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
                 inputX >> X2_vector[i](row, col);
            }
        }
    }

    // Close the input file
    inputX.close();

        // Create an array to store the matrices
    // Open the text file for reading
    std::ifstream inputV("Thermalised_VNis3.txt");
    if (!inputV.is_open()) {
        std::cerr << "Failed to open the V file." << std::endl;
        return 1;
    }


    // Read the values from the file and store them in the matrices
    for (int i = 0; i < dim; ++i) 
    {

        for (int row = 0; row < rows; ++row) 
        {
            for (int col = 0; col < cols; ++col) 
            {
                 inputV >> V2_vector[i](row, col);
            }
        }
    }

    // Close the input file
    inputV.close();

        // Create an array to store the matrices
    // Open the text file for reading
    std::ifstream inputA("Thermalised_ANis3.txt");
    if (!inputA.is_open()) {
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
                 inputA >> A2_vector[i](row, col);
            }
        }
    }

    // Close the input file
    inputA.close();

// If you want to test if the coords are loaded in right (They are)


    std::cout << "Start" << std::endl << "X";

    // Print out starting values
    for (matrix Matrix : X2_vector)
    {
        std::cout << Matrix << std::endl;
    }

    std::cout <<std::endl << "V";

    // Print out starting values
    for (matrix Matrix : V2_vector)
    {
        std::cout << Matrix << std::endl;
    }

    std::cout << std::endl << "A";

    // Print out starting values
    for (matrix Matrix : A2_vector)
    {
        std::cout << Matrix << std::endl;
    }



    matrix zero_matrix = matrix::Zero(rows, cols);

    matrix X2_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix V2_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix A2_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    
    // Perturb the system to get the two set of conditions, for 1 second before turning off the deformation.
    for (int j = 0; j < 1.0 / delta_t; ++j)
    {
        //for (int h = 0; h < dim; h++)
        //{
        //    std::cout << std::endl << X2_vector[h];
       // }
        // velocity Verlet 1 to get new positions from old positions, momentums and rate of change of momentums

        for (int i = 0; i < 9; ++i)
        {
            X2_vector_new[i] = X2_vector[i] + V2_vector[i] * delta_t + 0.5 * A2_vector[i] * delta_t * delta_t;
        }

        // Generate and store new A1, A2, A3, A4, A5, A6, A7, A8, and A9
        for (int i = 0; i < 9; ++i) 
        {
            A2_vector_new[i] = PerturbingAcceleration( i, X2_vector_new, rows, cols, g, c_1, c_2);
        }  
        
        // Use Velocity Verlet 2 to get new momentums
        for (int i = 0; i < 9; ++i)
        {
            V2_vector_new[i] = V2_vector[i] + 0.5 * (A2_vector_new[i] + A2_vector[i]) * delta_t;
        }
        
        // Copy elements from X_vector_new to X_vector
        std::memcpy(X2_vector, X2_vector_new, sizeof(X2_vector_new));  


        // Copy elements from X_vector_new to X_vector
        std::memcpy(V2_vector, V2_vector_new, sizeof(V2_vector_new));  

        // Copy elements from X_vector_new to X_vector
        std::memcpy(A2_vector, A2_vector_new, sizeof(A2_vector_new));  

        if (j % 1000 == 0)
        {
            //for (matrix el : V_vector)
            //{
            //    std::cout <<"Ideal " << el << std::endl;
            //}
            std::cout  << std::endl << gauss_law(X2_vector_new[0], X2_vector_new[1], X2_vector_new[2], X2_vector_new[3], X2_vector_new[4], X2_vector_new[5], X2_vector_new[6], X2_vector_new[7], X2_vector_new[8],
                                   V2_vector_new[0], V2_vector_new[1], V2_vector_new[2], V2_vector_new[3], V2_vector_new[4], V2_vector_new[5], V2_vector_new[6], V2_vector_new[7], V2_vector_new[8]);

            std::cout << std::endl;
            std::cout << "H" << std::setprecision(15) << H(1.0, 
                            X2_vector_new[0], X2_vector_new[1], X2_vector_new[2], X2_vector_new[3], X2_vector_new[4], X2_vector_new[5], X2_vector_new[6], X2_vector_new[7], X2_vector_new[8],
                          V2_vector_new[0], V2_vector_new[1], V2_vector_new[2], V2_vector_new[3], V2_vector_new[4], V2_vector_new[5], V2_vector_new[6], V2_vector_new[7], V2_vector_new[8]);


        }        

    }



    // Export initial X/V/A_vector to text files to be analysed in python.
    std:: fstream X2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D0-BranesNis3/perturbed_XNis3.txt", std:: ios:: out);
    X2_vector_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (matrix Matrix : X2_vector)
    {
        X2_vector_Export << Matrix << std::endl;
    }

    std:: fstream V2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D0-BranesNis3/perturbed_VNis3.txt", std:: ios:: out);
    V2_vector_Export << std::fixed << std::setprecision(15);
    for (matrix Matrix : V2_vector)
    {
        V2_vector_Export << Matrix << std::endl;
    }

    std:: fstream A2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D0-BranesNis3/perturbed_ANis3.txt", std:: ios:: out);
    A2_vector_Export << std::fixed << std::setprecision(15);
    // Print to text file
    for (matrix Matrix : A2_vector)
    {
        A2_vector_Export << Matrix << std::endl;
    }

    X2_vector_Export.close();
    V2_vector_Export.close();
    A2_vector_Export.close();


    return 0;
}