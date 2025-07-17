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
const double delta_t = 5e-4;
const double seconds_simulated = 200;

// Repeat simulation for 1000 seconds.
const int simulation_repetitions = seconds_simulated / delta_t;
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

// Cillian's Hamiltonian
double H(
    const double g,
    matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7, matrix X8, matrix X9,
    matrix V1, matrix V2, matrix V3, matrix V4, matrix V5, matrix V6, matrix V7, matrix V8, matrix V9)
{
    // Compute kinetic energy T
    R_or_C T = 1.0/(g * g * 2.0) * (V1*V1 + V2*V2 + V3*V3 + V4*V4 + V5*V5 + V6*V6 + V7*V7 + V8*V8 + V9*V9).trace();

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


    R_or_C U = - 1.0/(g * g *4.0) * commutator_sum.trace();

 
    return std::abs(T + U);
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

// Distance
double d_0(matrix* X1_vector, matrix* X2_vector)
{
    matrix delta_X = matrix::Zero(rows, cols);
    double sumX = 0.0;
    double argument = 0.0;

    for (int i = 0 ; i < dim ; i ++)
    {
        delta_X = X1_vector[i]-X2_vector[i];

        argument = (delta_X * delta_X).trace().real();

        sumX += argument;

    }

    return sqrt(sumX);
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


int main() 
{

    matrix zero_matrix = matrix::Zero(rows, cols);

    // Create  vectors to store the matrices
    matrix X1_vector[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix V1_vector[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix A1_vector[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix X2_vector[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix V2_vector[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix A2_vector[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};




    // Generate and store X1, X2, X3, X4, X5, X6, X7, X8, and X9
        // Create an array to store the matrices
    // Open the text file for reading
    std::ifstream inputX("thermalised_X.txt");
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
                 inputX >> X1_vector[i](row, col);
            }
        }
    }

    // Close the input file
    inputX.close();

        // Create an array to store the matrices
    // Open the text file for reading
    std::ifstream inputV("thermalised_V.txt");
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
    std::ifstream inputA("thermalised_A.txt");
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

 // Close the input file
    inputA.close();

    std::ifstream inputX2("perturbed_X.txt");
    if (!inputX2.is_open()) {
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
                 inputX2 >> X2_vector[i](row, col);
            }
        }
    }

    // Close the input file
    inputX2.close();

        // Create an array to store the matrices
    // Open the text file for reading
    std::ifstream inputV2("perturbed_V.txt");
    if (!inputV2.is_open()) {
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
                 inputV2 >> V2_vector[i](row, col);
            }
        }
    }

    // Close the input file
    inputV2.close();

        // Create an array to store the matrices
    // Open the text file for reading
    std::ifstream inputA2("perturbed_A.txt");
    if (!inputA2.is_open()) {
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
                 inputA2 >> A2_vector[i](row, col);
            }
        }
    }

    // Close the input file
    inputA2.close();

    matrix X1_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix X2_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix V1_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix V2_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix A1_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix A2_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix X2_vector_new_sprott[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};


    // Export initial X/V/A_vector to text files to be analysed in python.
    std:: fstream X1_vector_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/X_original_sprott_sim5.txt", std:: ios:: out);
    X1_vector_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (matrix Matrix : X1_vector)
    {
        X1_vector_Export << Matrix << std::endl;
    }

    std:: fstream V1_vector_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/V_original_sprott_sim5.txt", std:: ios:: out);
    V1_vector_Export << std::fixed << std::setprecision(15);
    for (matrix Matrix : V1_vector)
    {
        V1_vector_Export << Matrix << std::endl;
    }


    std:: fstream X2_vector_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/X_perturbed_1_sprott_sim5.txt", std:: ios:: out);
    X2_vector_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (matrix Matrix : X2_vector)
    {
        X2_vector_Export << Matrix << std::endl;
    }

    std:: fstream V2_vector_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/V_perturbed_1_sprott_sim5.txt", std:: ios:: out);
    V2_vector_Export << std::fixed << std::setprecision(15);
    for (matrix Matrix : V2_vector)
    {
        V2_vector_Export << Matrix << std::endl;
    }


    std::cout << X1_vector[0];




// Write simulation to advance the original coordinates by 1 second.


    for (int j = 0; j < 0.1 / delta_t; ++j)
    {

        // velocity Verlet 1 to get new positions from old positions, momentums and rate of change of momentums

        for (int i = 0; i < 9; ++i)
        {
            X1_vector_new[i] = X1_vector[i] + V1_vector[i] * delta_t + 0.5 * A1_vector[i] * delta_t * delta_t;
        }

        // Generate and store new A1, A2, A3, A4, A5, A6, A7, A8, and A9
        for (int i = 0; i < 9; ++i) 
        {
            A1_vector_new[i] = Acceleration( i, X1_vector_new, rows, cols, g);
        }  
        
        // Use Velocity Verlet 2 to get new momentums
        for (int i = 0; i < 9; ++i)
        {
            V1_vector_new[i] = V1_vector[i] + 0.5 * (A1_vector_new[i] + A1_vector[i]) * delta_t;

        }

        // Copy elements from X_vector_new to X_vector
        std::memcpy(X1_vector, X1_vector_new, sizeof(X1_vector_new));  


        // Copy elements from X_vector_new to X_vector
        std::memcpy(V1_vector, V1_vector_new, sizeof(V1_vector_new)); 

        // Copy elements from X_vector_new to X_vector
        std::memcpy(A1_vector, A1_vector_new, sizeof(A1_vector_new)); 

    }


    double distance_initial = d_0(X1_vector, X2_vector);

    double distance_new = 0.0;

    // Write simulation to simulate system
    for (int j = 0; j < seconds_simulated / delta_t; ++j)
    {

        // velocity Verlet 1 to get new positions from old positions, momentums and rate of change of momentums

        for (int i = 0; i < 9; ++i)
        {
            X1_vector_new[i] = X1_vector[i] + V1_vector[i] * delta_t + 0.5 * A1_vector[i] * delta_t * delta_t;


            X2_vector_new[i] = X2_vector[i] + V2_vector[i] * delta_t + 0.5 * A2_vector[i] * delta_t * delta_t;

        }
        // Generate and store new A1, A2, A3, A4, A5, A6, A7, A8, and A9
        for (int i = 0; i < 9; ++i) 
        {
            A1_vector_new[i] = Acceleration( i, X1_vector_new, rows, cols, g);

            A2_vector_new[i] = Acceleration( i, X2_vector_new, rows, cols, g);

        }  
        
        // Use Velocity Verlet 2 to get new momentums
        for (int i = 0; i < 9; ++i)
        {
            V1_vector_new[i] = V1_vector[i] + 0.5 * (A1_vector_new[i] + A1_vector[i]) * delta_t;

            V2_vector_new[i] = V2_vector[i] + 0.5 * (A2_vector_new[i] + A2_vector[i]) * delta_t;

        }
        

        distance_new = d_0(X1_vector_new, X2_vector_new);

        for (int i = 0; i < 9; i++)
        {
            X2_vector_new_sprott[i] = X1_vector_new[i] + distance_initial/distance_new * (X2_vector_new[i] - X1_vector_new[i]);
        }

        // Copy elements from X_vector_new to X_vector
        std::memcpy(X1_vector, X1_vector_new, sizeof(X1_vector_new));  

        // Copy elements from X_vector_new to X_vector
        std::memcpy(V1_vector, V1_vector_new, sizeof(V1_vector_new)); 

        // Copy elements from X_vector_new to X_vector
        std::memcpy(A1_vector, A1_vector_new, sizeof(A1_vector_new)); 


        // Copy elements from X_vector_new to X_vector
        std::memcpy(X2_vector, X2_vector_new_sprott, sizeof(X2_vector_new));  

        // Copy elements from X_vector_new to X_vector
        std::memcpy(V2_vector, V2_vector_new, sizeof(V2_vector_new));  

        // Copy elements from X_vector_new to X_vector
        std::memcpy(A2_vector, A2_vector_new, sizeof(A2_vector_new));  



        if (j % 10000 == 0)
        {
            std::cout << '\n' << "H" << std::setprecision(15) << H(g, 
                            X1_vector_new[0], X1_vector_new[1], X1_vector_new[2], X1_vector_new[3], X1_vector_new[4], X1_vector_new[5], X1_vector_new[6], X1_vector_new[7], X1_vector_new[8],
                          V1_vector_new[0], V1_vector_new[1], V1_vector_new[2], V1_vector_new[3], V1_vector_new[4], V1_vector_new[5], V1_vector_new[6], V1_vector_new[7], V1_vector_new[8]);

            std::cout << '\n' << j  <<"H2 :" << std::setprecision(15) << H(g, 
                            X2_vector_new[0], X2_vector_new[1], X2_vector_new[2], X2_vector_new[3], X2_vector_new[4], X2_vector_new[5], X2_vector_new[6], X2_vector_new[7], X2_vector_new[8],
                          V2_vector_new[0], V2_vector_new[1], V2_vector_new[2], V2_vector_new[3], V2_vector_new[4], V2_vector_new[5], V2_vector_new[6], V2_vector_new[7], V2_vector_new[8]);

        }

        if (j % 1000 == 0)
        {
            // gauss_law(X_vector_new[0], X_vector_new[1], X_vector_new[2], X_vector_new[3], X_vector_new[4], X_vector_new[5], X_vector_new[6], X_vector_new[7], X_vector_new[8],
            //                V_vector_new[0], V_vector_new[1], V_vector_new[2], V_vector_new[3], V_vector_new[4], V_vector_new[5], V_vector_new[6], V_vector_new[7], V_vector_new[8]);

            for ( matrix Matrix : X1_vector_new)
            {
                X1_vector_Export << Matrix << std::endl;
            }

            for ( matrix Matrix : V1_vector_new)
            {
                V1_vector_Export << Matrix << std::endl;
            }

            for ( matrix Matrix : X2_vector_new)
            {
                X2_vector_Export << Matrix << std::endl;
            }

            for ( matrix Matrix : V2_vector_new)
            {
                V2_vector_Export << Matrix << std::endl;
            }


        }

    }


    X2_vector_Export.close();
    V2_vector_Export.close();

    X1_vector_Export.close();
    V1_vector_Export.close();

    return 0;
} 