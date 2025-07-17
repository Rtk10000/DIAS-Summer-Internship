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

std::complex<long double> I = {0.0, 1.0};

long double c1 = 2.5;
long double c2 = 3.3;
long double c3 = 1.2;
long double c4 = -c1-c2-c3;

long double c12 = 1.1;
long double c13 = 2.2;
long double c23 = 2.0;
long double c14 = 0.7;
long double c24 = 3.1;
long double c34 = 0.9;

// Define timestep
const long double delta_t = 1e-3;
const long double seconds_thermalised = 10;

// Repeat simulation for 1000 seconds.
const int simulation_repetitions = seconds_thermalised / delta_t;
// Number of D0-Branes
const int N = 3;
const long double g = 1.0/sqrt(N);
const int rows = N;
const int cols = N;

// Dimension of space
const int dim = 9;

// l_or_d stands for long double or double, switch out double for long double if accuracy is low.
// R_or_C stands for Real or complex but means complex in this code, legacy from when I was doing real valued X elements in June.
typedef double l_or_d;
typedef std::complex<l_or_d> R_or_C;
typedef Eigen:: Matrix<std::complex<l_or_d>, N, N> matrix;
typedef Eigen:: Matrix<std::complex<l_or_d>, 1, N> row_vector;
typedef Eigen:: Matrix<std::complex<l_or_d>, N, 1> col_vector;

// Define the imaginary unit number i with the same symbol sympy uses.
R_or_C I(0,1);


matrix commutator(matrix A, matrix B)
{
    return A * B - B * A;
} 

// Initialising all Xs and Vs and As.
l_or_d X_scalars_vector[27] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
matrix X_matrices_vector[9] = {matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols)};
l_or_d V_scalars_vector[27] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
matrix V_matrices_vector[9] = {matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols)};
l_or_d A_scalars_vector[27] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
matrix A_matrices_vector[9] = {matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols)};


// Initialising all Zs and Z_dots.
// Z_scalars_vector = {Z12, Z13, Z21, Z23, Z31, Z32};
R_or_C C_zero(0.0, 0.0);
R_or_C Z_scalars_vector[6] = {C_zero, C_zero, C_zero, C_zero, C_zero, C_zero};
R_or_C Z_dot_scalars_vector[6] = {C_zero, C_zero, C_zero, C_zero, C_zero, C_zero};

// Z_row_vectors_vector = {Z14, Z24, Z34};
row_vector Z_row_vectors_vector[3] = {matrix::Zero(1, cols), matrix::Zero(1, cols), matrix::Zero(1, cols)};
row_vector Z_dot_row_vectors_vector[3] = {matrix::Zero(1, cols), matrix::Zero(1, cols), matrix::Zero(1, cols)};

// Z_row_vectors_vector = {Z41, Z42, Z43};
col_vector Z_col_vectors_vector[3] = {matrix::Zero(rows, 1), matrix::Zero(rows, 1), matrix::Zero(rows, 1)};
col_vector Z_dot_col_vectors_vector[3] = {matrix::Zero(rows, 1), matrix::Zero(rows, 1), matrix::Zero(rows, 1)};

// Initialising all phis and phi_dots


// Initialising all Hs 
matrix H(matrix Xa, matrix Xb, );

// Initialising all Gs


// Initialising all Fs


// Initialising all Ys

// Total energy
long double H( long double g, l_or_d* X_scalars_vector, matrix* X_matrices_vector, 
               l_or_d* V_scalars_vector, matrix* V_matrices_vector, l_or_d* A_scalars_vector,
               matrix* A_matrices_vector)
{
    l_or_d X_sum = 0.0;
    // Compute kinetic energy K
    for(l_or_d X: X_scalars_vector)
    {
        X_sum += X*X;
    }

    for(matrix X: X_matrices_vector)
    {
        X_sum += (X*X).trace();
    }

    l_or_d Z_sum = 0.0;

    for(R_or_C Z : Z_dot_scalars_vector)
    {
        Z_sum += Z*Z;
    }    
    
    for(col_vector Z : Z_dot_row_vectors_vector)
    {
        Z_sum += Z*Z.trace();
    }    
        
    for(col_vector Z : Z_dot_col_vectors_vector)
    {
        Z_sum += Z*Z.trace();
    }    
    l_or_d K = 0.5 * (  X_sum + Z_sum );

    l_or_d V_gauge = (  (H12_1_dag*H12_1) + (H12_2_dag*H12_2) + (H12_3_dag*H12_3) + 
                 (H21_1_dag*H21_1) + (H21_2_dag*H21_2) + (H21_3_dag*H21_3) + 
                 (H13_1_dag*H13_1) + (H13_2_dag*H13_2) + (H13_3_dag*H13_3) + 
                 (H31_1_dag*H31_1) + (H31_2_dag*H31_2) + (H31_3_dag*H31_3) + 
                 (H23_1_dag*H23_1) + (H23_2_dag*H23_2) + (H23_3_dag*H23_3) + 
                 (H32_1_dag*H32_1) + (H32_2_dag*H32_2) + (H32_3_dag*H32_3) + 
                 ((H14_1_dag*H14_1)).trace() + ((H14_2_dag*H14_2)).trace() + ((H14_3_dag*H14_3)).trace()+ 
                 ((H41_1_dag*H41_1)).trace() + ((H41_2_dag*H41_2)).trace() + ((H41_3_dag*H41_3)).trace()+ 
                 ((H24_1_dag*H24_1)).trace() + ((H24_2_dag*H24_2)).trace() + ((H24_3_dag*H24_3)).trace()+ 
                 ((H42_1_dag*H42_1)).trace() + ((H42_2_dag*H42_2)).trace() + ((H42_3_dag*H42_3)).trace()+ 
                 ((H34_1_dag*H34_1)).trace() + ((H34_2_dag*H34_2)).trace() + ((H34_3_dag*H34_3)).trace()+ 
                 ((H43_1_dag*H43_1)).trace() + ((H43_2_dag*H43_2)).trace() + ((H43_3_dag*H43_3)).trace()+
                 # \sum( \sum( \sum( Tr[ Dagger(Y^k_aA) * Y^k_aA] ) ) ) term
                 (Y1_11_dag * Y1_11) + (Y1_11_dag * Y1_11) + (Y1_11_dag * Y1_11) + 
                 (Y1_21_dag * Y1_21) + (Y1_22_dag * Y1_22) + (Y1_23_dag * Y1_23) +
                 (Y1_31_dag * Y1_31) + (Y1_32_dag * Y1_32) + (Y1_32_dag * Y1_32) + 
                 (Y2_11_dag * Y2_11) + (Y2_11_dag * Y2_11) + (Y2_11_dag * Y2_11) + 
                 (Y2_21_dag * Y2_21) + (Y2_22_dag * Y2_22) + (Y2_23_dag * Y2_23) +
                 (Y2_31_dag * Y2_31) + (Y2_32_dag * Y2_32) + (Y2_32_dag * Y2_32) + 
                 (Y3_11_dag * Y3_11) + (Y3_11_dag * Y3_11) + (Y3_11_dag * Y3_11) + 
                 (Y3_21_dag * Y3_21) + (Y3_22_dag * Y3_22) + (Y3_23_dag * Y3_23) +
                 (Y3_31_dag * Y3_31) + (Y3_32_dag * Y3_32) + (Y3_32_dag * Y3_32) +
                 ((Y4_11_dag * Y4_11)).trace() + ((Y4_11_dag * Y4_11)).trace() + (Y4_11_dag * Y4_11).trace() + 
                 ((Y4_21_dag * Y4_21)).trace() + ((Y4_22_dag * Y4_22)).trace() + (Y4_23_dag * Y4_23).trace() +
                 ((Y4_31_dag * Y4_31)).trace() + ((Y4_32_dag * Y4_32)).trace() + (Y4_32_dag * Y4_32).trace() +
                 1/4 *( X1_12*X1_12 + X1_21*X1_21 + X1_13*X1_13 + X1_31*X1_31 +
                        X1_23*X1_23 + X1_32*X1_32 +
                        X2_12*X2_12 + X2_21*X2_21 + X2_13*X2_13 + X2_31*X2_31 +
                        X2_23*X2_23 + X2_32*X2_32 +
                        X3_12*X3_12 + X3_21*X3_21 + X3_13*X3_13 + X3_31*X3_31 +
                        X3_23*X3_23 + X3_32*X3_32 +
                        ( X4_12*X4_12 ).trace() + ( X4_21*X4_21 ).trace() + ( X4_13*X4_13 ).trace() + 
                        ( X4_31*X4_31 ).trace() + ( X4_23*X4_23 ).trace() + ( X4_32*X4_32 ).trace()
                      )
               );

    l_or_d V_D_arg_k1 = (Z12 * Z12_dag + Z13 * Z13_dag + Trace(Z14 * Z14_dag) -
              ( Z21_dag * Z21 + Z31_dag * Z31 + Trace(Z41_dag * Z41) ) + 
              ( (phi1_1*phi1_1_dag).expand() - (phi1_1_dag*phi1_1).expand() + (phi1_2*phi1_2_dag).expand() - 
                  (phi1_2_dag*phi1_2).expand() + (phi1_3*phi1_3_dag).expand() - (phi1_3_dag*phi1_3).expand() ) - 
              c1/(g*g) );

    l_or_d V_D_arg_k2 = (Z21 * Z21_dag + Z23 * Z23_dag + Trace(Z24 * Z24_dag) - 
              ( Z12_dag * Z12 + Z32_dag * Z32 + Trace(Z42_dag * Z42) ) + 
              ( (phi2_1*phi2_1_dag).expand() - (phi2_1_dag*phi2_1).expand() + (phi2_2*phi2_2_dag).expand() - 
                (phi2_2_dag*phi2_2).expand() + (phi2_3*phi2_3_dag).expand() - (phi2_3_dag*phi2_3).expand() ) - 
              c2/(g*g) );

    l_or_d V_D_arg_k3 = (Z31 * Z31_dag + Z32 * Z32_dag + Trace(Z34 * Z34_dag) - 
              ( Z13_dag * Z13 + Z23_dag * Z23 + Trace(Z43_dag * Z43) ) + 
              ( (phi3_1*phi3_1_dag).expand() - (phi3_1_dag*phi3_1).expand() + (phi3_2*phi3_2_dag).expand() - 
                (phi3_2_dag*phi3_2).expand() + (phi3_3*phi3_3_dag).expand() - (phi3_3_dag*phi3_3).expand() ) - 
              c3/(g*g) );

    l_or_d V_D_arg_k4 = (Z41 * Z41_dag + Z42 * Z42_dag + Z43 * Z43_dag + 
              ( Z14_dag * Z14 + Z24_dag * Z24 + Z34_dag * Z34 ) + 
              ( (phi4_1*phi4_1_dag).expand() - (phi4_1_dag*phi4_1).expand() + (phi4_2*phi4_2_dag).expand() - 
                (phi4_2_dag*phi4_2).expand() + (phi4_3*phi4_3_dag).expand() - (phi4_3_dag*phi4_3).expand() ) - 
              c4/(g*g)*Identity(N) );


    l_or_d V_D = 1/2 * ((V_D_arg_k1 * V_D_arg_k1).expand() + (V_D_arg_k2 * V_D_arg_k2).expand() + 
             (V_D_arg_k3 * V_D_arg_k3).expand() + Trace(V_D_arg_k4 * V_D_arg_k4).expand() );

    l_or_d V_F = ( Abs(F12)*Abs(F12) + Abs(F21)*Abs(F21) + Abs(F13)*Abs(F13) + Abs(F31)*Abs(F31) + Abs(F32)*Abs(F32) + Abs(F23)*Abs(F23) +
                   F14 * F14_dag + F24 * F24_dag + F34 * F34_dag + 
                   Trace( F41_dag * F41 ) + Trace( F42_dag * F42) + Trace( F43_dag * F43 )+
                   Abs(G12)*Abs(G12) + Abs(G21)*Abs(G21) + Abs(G13)*Abs(G13) + Abs(G31)*Abs(G31) + Abs(G32)*Abs(G32) + Abs(G23)*Abs(G23) + 
                   Trace( (G41_dag * G41) ).expand() + Trace( (G14 * G14_dag) ).expand() + Trace( (G42_dag * G42) ).expand() + 
                   Trace( (G24 * G24_dag) ).expand() + Trace( (G43_dag * G43) ).expand() + Trace( (G34 * G34_dag) ).expand() )

    return K - 0.5*g*g*(V_gauge, V_D, V_F);
}

// Cillian's Gauss' law
matrix gauss_law(
    matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7, matrix X8, matrix X9,
    matrix V1, matrix V2, matrix V3, matrix V4, matrix V5, matrix V6, matrix V7, matrix V8, matrix V9)
{
    matrix result =  matrix::Zero(rows, cols);
    for (int i = 0; i < 9; i ++)
    {
        result = commutator(X1,V1) + commutator(X2,V2) + commutator(X3,V3) + commutator(X4,V4) + commutator(X5,V5) +
        commutator(X6,V6) + commutator(X7,V7) + commutator(X8,V8) + commutator(X9,V9);
    }
    return result;
}


// Random Hermitian matrix

// Acceleration of each coordinate matrix
matrix Acceleration(const int j, matrix* X_vector, int rows, int cols, const long double g)
{

    matrix temp_commutator = matrix::Zero(rows, cols);
    matrix commutator_sum = matrix::Zero(rows, cols);

    matrix X = X_vector[j];
    for (int i = 0; i < dim; ++i)
    {
        if (i != j)
        {  
            temp_commutator = commutator(X_vector[i], commutator(X, X_vector[i]));
            

            commutator_sum += temp_commutator;
            
        }   
    // Comment out the line below to go back to the lagrangian e.q.m. With no -1/(g^2)
    // commutator_sum = -1.0 / (g * g) * commutator_sum;
    }

    return g*g*commutator_sum;
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
    std::ifstream inputX("initial_X.txt");
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
                X_vector[i](row, col) = X_vector[i](row, col)/g;
            }
        }
    }

    // Close the input file
    inputX.close();

    

    // Generate and store A1, A2, A3, A4, A5, A6, A7, A8, and A9
    for (int i = 0; i < dim; ++i) 
    {
        A_vector[i] = Acceleration(i, X_vector, rows, cols, g);
    }  


    // Create  vectors to store the new matrices
    matrix X_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix V_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix A_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};

    // Defining all Xs and X_dots
    R_or_C X1_1(X_scalars_vector[0], 0);
    R_or_C X1_2(X_scalars_vector[1], 0);
    R_or_C X1_3(X_scalars_vector[2], 0);
    R_or_C X1_4(X_scalars_vector[3], 0);
    R_or_C X1_5(X_scalars_vector[4], 0);
    R_or_C X1_6(X_scalars_vector[5], 0);
    R_or_C X1_7(X_scalars_vector[6], 0);
    R_or_C X1_8(X_scalars_vector[7], 0);
    R_or_C X1_9(X_scalars_vector[8], 0);

    R_or_C X2_1(X_scalars_vector[9], 0);
    R_or_C X2_2(X_scalars_vector[10], 0);
    R_or_C X2_3(X_scalars_vector[11], 0);
    R_or_C X2_4(X_scalars_vector[12], 0);
    R_or_C X2_5(X_scalars_vector[13], 0);
    R_or_C X2_6(X_scalars_vector[14], 0);
    R_or_C X2_7(X_scalars_vector[15], 0);
    R_or_C X2_8(X_scalars_vector[16], 0);
    R_or_C X2_9(X_scalars_vector[17], 0);

    R_or_C X3_1(X_scalars_vector[18], 0);
    R_or_C X3_2(X_scalars_vector[19], 0);
    R_or_C X3_3(X_scalars_vector[20], 0);
    R_or_C X3_4(X_scalars_vector[21], 0);
    R_or_C X3_5(X_scalars_vector[22], 0);
    R_or_C X3_6(X_scalars_vector[23], 0);
    R_or_C X3_7(X_scalars_vector[24], 0);
    R_or_C X3_8(X_scalars_vector[25], 0);
    R_or_C X3_9(X_scalars_vector[26], 0);

    matrix X4_1(X_matrices_vector[0], 0);
    matrix X4_2(X_matrices_vector[1], 0);
    matrix X4_3(X_matrices_vector[2], 0);
    matrix X4_4(X_matrices_vector[3], 0);
    matrix X4_5(X_matrices_vector[4], 0);
    matrix X4_6(X_matrices_vector[5], 0);
    matrix X4_7(X_matrices_vector[6], 0);
    matrix X4_8(X_matrices_vector[7], 0);
    matrix X4_9(X_matrices_vector[8], 0);


    // Defining all phis and phi_dots
    R_or_C coeff(1.0/sqrt(2.0), 0);
    R_or_C phi1_1 = coeff*(X1_4+I*X1_5);
    R_or_C phi1_2 = coeff*(X1_6+I*X1_7);
    R_or_C phi1_3 = coeff*(X1_8+I*X1_9);

    R_or_C phi1_1_dag = coeff*(X1_4-I*X1_5);
    R_or_C phi1_2_dag = coeff*(X1_6-I*X1_7);
    R_or_C phi1_3_dag = coeff*(X1_8-I*X1_9);

    R_or_C phi2_1 = coeff*(X2_4+I*X2_5);
    R_or_C phi2_2 = coeff*(X2_6+I*X2_7);
    R_or_C phi2_3 = coeff*(X2_8+I*X2_9);

    R_or_C phi2_1_dag = coeff*(X2_4-I*X2_5);
    R_or_C phi2_2_dag = coeff*(X2_6-I*X2_7);
    R_or_C phi2_3_dag = coeff*(X2_8-I*X2_9);

    R_or_C phi3_1 = coeff*(X3_4+I*X3_5);
    R_or_C phi3_2 = coeff*(X3_6+I*X3_7);
    R_or_C phi3_3 = coeff*(X3_8+I*X3_9);

    R_or_C phi3_1_dag = coeff*(X3_4-I*X3_5);
    R_or_C phi3_2_dag = coeff*(X3_6-I*X3_7);
    R_or_C phi3_3_dag = coeff*(X3_8-I*X3_9);

    matrix phi4_1 = coeff*(X4_4+I*X4_5);
    matrix phi4_2 = coeff*(X4_6+I*X4_7);
    matrix phi4_3 = coeff*(X4_8+I*X4_9);

    matrix phi4_1_dag = coeff*(X4_4-I*X4_5);
    matrix phi4_2_dag = coeff*(X4_6-I*X4_7);
    matrix phi4_3_dag = coeff*(X4_8-I*X4_9);

    // Defining phi12 and phi12_dag
    R_or_C phi12 = phi1_3 - phi2_3;
    R_or_C phi23 = phi2_1 - phi3_1;
    R_or_C phi31 = phi3_2 - phi1_2;

    R_or_C phi12_dag = phi1_3_dag - phi2_3_dag;
    R_or_C phi23_dag = phi2_1_dag - phi3_1_dag;
    R_or_C phi31_dag = phi3_2_dag - phi1_2_dag;

    // Defining Ys and Y_dots and _dags of them
R_or_C Y1_11 = X1_1 * phi1_1 - phi1_1  * X1_1;
R_or_C Y1_12 = X1_1 * phi1_2 - phi1_2  * X1_1;
R_or_C Y1_13 = X1_1 * phi1_3 - phi1_3  * X1_1;
R_or_C Y1_21 = X1_2 * phi1_1 - phi1_1  * X1_2;
R_or_C Y1_22 = X1_2 * phi1_2 - phi1_2  * X1_2;
R_or_C Y1_23 = X1_2 * phi1_3 - phi1_3  * X1_2;
R_or_C Y1_31 = X1_3 * phi1_1 - phi1_1  * X1_3;
R_or_C Y1_32 = X1_3 * phi1_2 - phi1_2  * X1_3;
R_or_C Y1_33 = X1_3 * phi1_3 - phi1_3  * X1_3;

R_or_C Y2_11 = X2_1 * phi2_1 - phi2_1  * X2_1;
R_or_C Y2_12 = X2_1 * phi2_2 - phi2_2  * X2_1;
R_or_C Y2_13 = X2_1 * phi2_3 - phi2_3  * X2_1;
R_or_C Y2_21 = X2_2 * phi2_1 - phi2_1  * X2_2;
R_or_C Y2_22 = X2_2 * phi2_2 - phi2_2  * X2_2;
R_or_C Y2_23 = X2_2 * phi2_3 - phi2_3  * X2_2;
R_or_C Y2_31 = X2_3 * phi2_1 - phi2_1  * X2_3;
R_or_C Y2_32 = X2_3 * phi2_2 - phi2_2  * X2_3;
R_or_C Y2_33 = X2_3 * phi2_3 - phi2_3  * X2_3;

R_or_C Y3_11 = X3_1 * phi3_1 - phi3_1  * X3_1 ;
R_or_C Y3_12 = X3_1 * phi3_2 - phi3_2  * X3_1 ;
R_or_C Y3_13 = X3_1 * phi3_3 - phi3_3  * X3_1 ;
R_or_C Y3_21 = X3_2 * phi3_1 - phi3_1  * X3_2 ;
R_or_C Y3_22 = X3_2 * phi3_2 - phi3_2  * X3_2 ;
R_or_C Y3_23 = X3_2 * phi3_3 - phi3_3  * X3_2 ;
R_or_C Y3_31 = X3_3 * phi3_1 - phi3_1  * X3_3 ;
R_or_C Y3_32 = X3_3 * phi3_2 - phi3_2  * X3_3 ;
R_or_C Y3_33 = X3_3 * phi3_3 - phi3_3  * X3_3 ;

matrix Y4_11 = ( X4_1 * phi4_1) - ( phi4_1   * X4_1);
matrix Y4_12 = ( X4_1 * phi4_2) - ( phi4_2   * X4_1);
matrix Y4_13 = ( X4_1 * phi4_3) - ( phi4_3   * X4_1);
matrix Y4_21 = ( X4_2 * phi4_1) - ( phi4_1   * X4_2);
matrix Y4_22 = ( X4_2 * phi4_2) - ( phi4_2   * X4_2);
matrix Y4_23 = ( X4_2 * phi4_3) - ( phi4_3   * X4_2);
matrix Y4_31 = ( X4_3 * phi4_1) - ( phi4_1   * X4_3);
matrix Y4_32 = ( X4_3 * phi4_2) - ( phi4_2   * X4_3);
matrix Y4_33 = ( X4_3 * phi4_3) - ( phi4_3   * X4_3);

R_or_C Y1_11_dag = X1_1 * phi1_1_dag - phi1_1_dag * X1_1;
R_or_C Y1_12_dag = X1_1 * phi1_2_dag - phi1_2_dag * X1_1;
R_or_C Y1_13_dag = X1_1 * phi1_3_dag - phi1_3_dag * X1_1;
R_or_C Y1_21_dag = X1_2 * phi1_1_dag - phi1_1_dag * X1_2;
R_or_C Y1_22_dag = X1_2 * phi1_2_dag - phi1_2_dag * X1_2;
R_or_C Y1_23_dag = X1_2 * phi1_3_dag - phi1_3_dag * X1_2;
R_or_C Y1_31_dag = X1_3 * phi1_1_dag - phi1_1_dag * X1_3;
R_or_C Y1_32_dag = X1_3 * phi1_2_dag - phi1_2_dag * X1_3;
R_or_C Y1_33_dag = X1_3 * phi1_3_dag - phi1_3_dag * X1_3;

R_or_C Y2_11_dag = X2_1 * phi2_1_dag - phi2_1_dag * X2_1;
R_or_C Y2_12_dag = X2_1 * phi2_2_dag - phi2_2_dag * X2_1;
R_or_C Y2_13_dag = X2_1 * phi2_3_dag - phi2_3_dag * X2_1;
R_or_C Y2_21_dag = X2_2 * phi2_1_dag - phi2_1_dag * X2_2;
R_or_C Y2_22_dag = X2_2 * phi2_2_dag - phi2_2_dag * X2_2;
R_or_C Y2_23_dag = X2_2 * phi2_3_dag - phi2_3_dag * X2_2;
R_or_C Y2_31_dag = X2_3 * phi2_1_dag - phi2_1_dag * X2_3;
R_or_C Y2_32_dag = X2_3 * phi2_2_dag - phi2_2_dag * X2_3;
R_or_C Y2_33_dag = X2_3 * phi2_3_dag - phi2_3_dag * X2_3;

R_or_C Y3_11_dag = X3_1 * phi3_1_dag - phi3_1_dag * X3_1;
R_or_C Y3_12_dag = X3_1 * phi3_2_dag - phi3_2_dag * X3_1;
R_or_C Y3_13_dag = X3_1 * phi3_3_dag - phi3_3_dag * X3_1;
R_or_C Y3_21_dag = X3_2 * phi3_1_dag - phi3_1_dag * X3_2;
R_or_C Y3_22_dag = X3_2 * phi3_2_dag - phi3_2_dag * X3_2;
R_or_C Y3_23_dag = X3_2 * phi3_3_dag - phi3_3_dag * X3_2;
R_or_C Y3_31_dag = X3_3 * phi3_1_dag - phi3_1_dag * X3_3;
R_or_C Y3_32_dag = X3_3 * phi3_2_dag - phi3_2_dag * X3_3;
R_or_C Y3_33_dag = X3_3 * phi3_3_dag - phi3_3_dag * X3_3;

matrix Y4_11_dag = ( X4_1 * phi4_1_dag) - (phi4_1_dag * X4_1);
matrix Y4_12_dag = ( X4_1 * phi4_2_dag) - (phi4_2_dag * X4_1);
matrix Y4_13_dag = ( X4_1 * phi4_3_dag) - (phi4_3_dag * X4_1);
matrix Y4_21_dag = ( X4_2 * phi4_1_dag) - (phi4_1_dag * X4_2);
matrix Y4_22_dag = ( X4_2 * phi4_2_dag) - (phi4_2_dag * X4_2);
matrix Y4_23_dag = ( X4_2 * phi4_3_dag) - (phi4_3_dag * X4_2);
matrix Y4_31_dag = ( X4_3 * phi4_1_dag) - (phi4_1_dag * X4_3);
matrix Y4_32_dag = ( X4_3 * phi4_2_dag) - (phi4_2_dag * X4_3);
matrix Y4_33_dag = ( X4_3 * phi4_3_dag) - (phi4_3_dag * X4_3);

    // Defining Fs and F_dags
    R_or_C F12 = Z12*Z21 + c12/(g*g);
    R_or_C F21 = Z21*Z12 + c12/(g*g);
    R_or_C F23 = Z23*Z32 + c23/(g*g);
    R_or_C F32 = Z32*Z23 + c23/(g*g);
    R_or_C F13 = Z13*Z31 + c13/(g*g);
    R_or_C F31 = Z31*Z13 + c13/(g*g);

    matrix F14 = Trace( MatMul(Z14,Z41) ) + c14*N/(g*g);
    matrix F41 = MatMul(Z41,Z14) + c14*Identity(N)/(g*g) + (phi4_2*phi4_3 - phi4_3*phi4_2);
    matrix F24 = Trace( MatMul(Z24,Z42) ) + c24*N/(g*g);
    matrix F42 = MatMul(Z42,Z24) + c24*Identity(N)/(g*g) + (phi4_3*phi4_1 - phi4_1*phi4_3);
    matrix F34 = Trace( MatMul(Z34,Z43) ) + c34*N/(g*g);
    matrix F43 = MatMul(Z43,Z34) + c34*Identity(N)/(g*g) - (phi4_2*phi4_1 - phi4_1*phi4_2);

    R_or_C F12_dag = Z12_dag*Z21_dag + c12/(g*g);
    R_or_C F21_dag = Z21_dag*Z12_dag + c12/(g*g);
    R_or_C F23_dag = Z23_dag*Z32_dag + c23/(g*g);
    R_or_C F32_dag = Z32_dag*Z23_dag + c23/(g*g);
    R_or_C F13_dag = Z13_dag*Z31_dag + c13/(g*g);
    R_or_C F31_dag = Z31_dag*Z13_dag + c13/(g*g);

    R_or_C F14_dag = Dagger(Trace( MatMul(Z14,Z41) ) )+ c14*N/(g*g);
    matrix F41_dag = Dagger(MatMul(Z41,Z14) )+ c14*Identity(N)/(g*g) + (-phi4_2_dag*phi4_3_dag + phi4_3_dag*phi4_2_dag);
    R_or_C F24_dag = Dagger(Trace( MatMul(Z24,Z42) ) )+ c24*N/(g*g);
    matrix F42_dag = Dagger(MatMul(Z42,Z24)) + c24*Identity(N)/(g*g) + (-phi4_3_dag*phi4_1_dag + phi4_1_dag*phi4_3_dag);
    R_or_C F34_dag = Dagger(Trace( MatMul(Z34,Z43) )) + c34*N/(g*g);
    matrix F43_dag = Dagger(MatMul(Z43,Z34) )+ c34*Identity(N)/(g*g) - (-phi4_2_dag*phi4_1_dag + phi4_1_dag*phi4_2_dag);

    // Defining Gs and G_dags
    R_or_C G21 = phi12 * Z21 + Z23*Z31 + Trace( MatMul(Z24,Z41) );
    R_or_C G12 = phi12 * Z12 + Z13*Z32 + Trace( MatMul(Z14,Z42) );
     
    R_or_C G31 = phi31 * Z31 + Z32*Z21 - Trace( MatMul(Z34,Z41) );
    R_or_C G13 = phi31 * Z13 + Z12*Z23 + Trace( MatMul(Z14,Z43) );
 
    R_or_C G32 = phi23 * Z32 + Z31*Z12 + Trace( MatMul(Z34,Z42) );
    R_or_C G23 = phi23 * Z23 + Z21*Z13 + Trace( MatMul(Z24,Z43) );

    matrix G41 = -phi4_1 * Z41 + Z42*Z21 + Z43*Z31;
    matrix G14 = -Z14 * phi4_1 + Z12*Z24 - Z13*Z34;
    
    matrix G42 = -phi4_2 * Z42 + Z43*Z32 + Z41*Z12;
    matrix G24 = -Z24 * phi4_2 + Z21*Z14 + Z23*Z34;

    matrix G43 = -phi4_3 * Z43 - Z41*Z13 + Z42*Z23;
    matrix G34 = -Z34 * phi4_3 + Z31*Z14 + Z32*Z24;

    R_or_C G21_dag = phi12_dag * Z21_dag + Z23_dag*Z31_dag + Dagger(Trace( MatMul(Z24,Z41) ));
    R_or_C G12_dag = phi12_dag * Z12_dag + Z13_dag*Z32_dag + Dagger(Trace( MatMul(Z14,Z42) ));

    R_or_C G31_dag = phi31_dag * Z31_dag + Z32_dag*Z21_dag - Dagger(Trace( MatMul(Z34,Z41) ));
    R_or_C G13_dag = phi31_dag * Z13_dag + Z12_dag*Z23_dag + Dagger(Trace( MatMul(Z14,Z43) ));

    R_or_C G32_dag = phi23_dag * Z32_dag + Z31_dag*Z12_dag + Dagger(Trace( MatMul(Z34,Z42) ));
    R_or_C G23_dag = phi23_dag * Z23_dag + Z21_dag*Z13_dag + Dagger(Trace( MatMul(Z24,Z43) ));

    row_vector G41_dag = -Z41_dag * phi4_1_dag + Z42_dag*Z21_dag + Z43_dag*Z31_dag;
    col_vector G14_dag = -phi4_1_dag * Z14_dag + Z12_dag*Z24_dag - Z13_dag*Z34_dag;

    row_vector G42_dag = -Z42_dag * phi4_2_dag + Z43_dag*Z32_dag + Z41_dag*Z12_dag;
    col_vector G24_dag = -phi4_2_dag * Z24_dag + Z21_dag*Z14_dag + Z23_dag*Z34_dag;

    row_vector G43_dag = -Z43_dag * phi4_3_dag - Z41_dag*Z13_dag + Z42_dag*Z23_dag;
    col_vector G34_dag = -phi4_3_dag * Z34_dag + Z31_dag*Z14_dag + Z32_dag*Z24_dag;

    // Defining Hs and H_dags
    R_or_C H12_1 = X1_1*Z12 - Z12*X2_1;
    R_or_C H12_2 = X1_2*Z12 - Z12*X2_2;
    R_or_C H12_3 = X1_3*Z12 - Z12*X2_3;

    R_or_C H21_1 = X2_1*Z21 - Z21*X1_1;
    R_or_C H21_2 = X2_2*Z21 - Z21*X1_2;
    R_or_C H21_3 = X2_3*Z21 - Z21*X1_3;

    R_or_C H13_1 = X1_1*Z13 - Z13*X3_1;
    R_or_C H13_2 = X1_2*Z13 - Z13*X3_2;
    R_or_C H13_3 = X1_3*Z13 - Z13*X3_3;

    R_or_C H31_1 = X3_1*Z31 - Z31*X1_1;
    R_or_C H31_2 = X3_2*Z31 - Z31*X1_2;
    R_or_C H31_3 = X3_3*Z31 - Z31*X1_3;

    col_vector H14_1 = X1_1*Z14 - Z14*X4_1;
    col_vector H14_2 = X1_2*Z14 - Z14*X4_2;
    col_vector H14_3 = X1_3*Z14 - Z14*X4_3;

    row_vector H41_1 = X4_1*Z41 - Z41*X1_1;
    row_vector H41_2 = X4_2*Z41 - Z41*X1_2;
    row_vector H41_3 = X4_3*Z41 - Z41*X1_3;

    R_or_C H32_1 = X3_1*Z32 - Z32*X2_1;
    R_or_C H32_2 = X3_2*Z32 - Z32*X2_2;
    R_or_C H32_3 = X3_3*Z32 - Z32*X2_3;

    R_or_C H23_1 = X2_1*Z23 - Z23*X3_1;
    R_or_C H23_2 = X2_2*Z23 - Z23*X3_2;
    R_or_C H23_3 = X2_3*Z23 - Z23*X3_3;

    row_vector H34_1 = X3_1*Z34 - Z34*X4_1;
    row_vector H34_2 = X3_2*Z34 - Z34*X4_2;
    row_vector H34_3 = X3_3*Z34 - Z34*X4_3;

    col_vector H43_1 = X4_1*Z43 - Z43*X3_1;
    col_vector H43_2 = X4_2*Z43 - Z43*X3_2;
    col_vector H43_3 = X4_3*Z43 - Z43*X3_3;

    row_vector H24_1 = X2_1*Z24 - Z24*X4_1;
    row_vector H24_2 = X2_2*Z24 - Z24*X4_2;
    row_vector H24_3 = X2_3*Z24 - Z24*X4_3;

    col_vector H42_1 = X4_1*Z42 - Z42*X2_1;
    col_vector H42_2 = X4_2*Z42 - Z42*X2_2;
    col_vector H42_3 = X4_3*Z42 - Z42*X2_3;

    R_or_C H12_1_dag = X1_1*Z12_dag - Z12_dag*X2_1;
    R_or_C H12_2_dag = X1_2*Z12_dag - Z12_dag*X2_2;
    R_or_C H12_3_dag = X1_3*Z12_dag - Z12_dag*X2_3;

    R_or_C H21_1_dag = X2_1*Z21_dag - Z21_dag*X1_1;
    R_or_C H21_2_dag = X2_2*Z21_dag - Z21_dag*X1_2;
    R_or_C H21_3_dag = X2_3*Z21_dag - Z21_dag*X1_3;

    R_or_C H13_1_dag = X1_1*Z13_dag - Z13_dag*X3_1;
    R_or_C H13_2_dag = X1_2*Z13_dag - Z13_dag*X3_2;
    R_or_C H13_3_dag = X1_3*Z13_dag - Z13_dag*X3_3;

    R_or_C H31_1_dag = X3_1*Z31_dag - Z31_dag*X1_1;
    R_or_C H31_2_dag = X3_2*Z31_dag - Z31_dag*X1_2;
    R_or_C H31_3_dag = X3_3*Z31_dag - Z31_dag*X1_3;

    row_vector H14_1_dag = Z14_dag*X1_1 - X4_1*Z14_dag;
    row_vector H14_2_dag = Z14_dag*X1_2 - X4_2*Z14_dag;
    row_vector H14_3_dag = Z14_dag*X1_3 - X4_3*Z14_dag;

    col_vector H41_1_dag = Z41_dag*X4_1 - X1_1*Z41_dag;
    col_vector H41_2_dag = Z41_dag*X4_2 - X1_2*Z41_dag;
    col_vector H41_3_dag = Z41_dag*X4_3 - X1_3*Z41_dag;

    R_or_C H32_1_dag = Z32_dag*X3_1 - X2_1*Z32_dag;
    R_or_C H32_2_dag = Z32_dag*X3_2 - X2_2*Z32_dag;
    R_or_C H32_3_dag = Z32_dag*X3_3 - X2_3*Z32_dag;

    R_or_C H23_1_dag = Z23_dag*X2_1 - X3_1*Z23_dag;
    R_or_C H23_2_dag = Z23_dag*X2_2 - X3_2*Z23_dag;
    R_or_C H23_3_dag = Z23_dag*X2_3 - X3_3*Z23_dag;

    row_vector H34_1_dag = Z34_dag*X3_1 - X4_1*Z34_dag;
    row_vector H34_2_dag = Z34_dag*X3_2 - X4_2*Z34_dag;
    row_vector H34_3_dag = Z34_dag*X3_3 - X4_3*Z34_dag;

    col_vector H43_1_dag = Z43_dag*X4_1 - X3_1*Z43_dag;
    col_vector H43_2_dag = Z43_dag*X4_2 - X3_2*Z43_dag;
    col_vector H43_3_dag = Z43_dag*X4_3 - X3_3*Z43_dag;

    row_vector H24_1_dag = Z24_dag*X2_1 - X4_1*Z24_dag;
    row_vector H24_2_dag = Z24_dag*X2_2 - X4_2*Z24_dag;
    row_vector H24_3_dag = Z24_dag*X2_3 - X4_3*Z24_dag;

    col_vector H42_1_dag = Z42_dag*X4_1 - X2_1*Z42_dag;
    col_vector H42_2_dag = Z42_dag*X4_2 - X2_2*Z42_dag;
    col_vector H42_3_dag = Z42_dag*X4_3 - X2_3*Z42_dag;

    // Xk_ab 
    R_or_C X1_12 = X1_1*X1_2 - X1_2*X1_1;
    R_or_C X1_21 = X1_2*X1_1 - X1_1*X1_2;
    R_or_C X1_13 = X1_1*X1_3 - X1_3*X1_1;
    R_or_C X1_31 = X1_3*X1_1 - X1_1*X1_3;
    R_or_C X1_23 = X1_2*X1_3 - X1_3*X1_2;
    R_or_C X1_32 = X1_3*X1_2 - X1_2*X1_3;

    R_or_C X2_12 = X2_1*X2_2 - X2_2*X2_1;
    R_or_C X2_21 = X2_2*X2_1 - X2_1*X2_2;
    R_or_C X2_13 = X2_1*X2_3 - X2_3*X2_1;
    R_or_C X2_31 = X2_3*X2_1 - X2_1*X2_3;
    R_or_C X2_23 = X2_2*X2_3 - X2_3*X2_2;
    R_or_C X2_32 = X2_3*X2_2 - X2_2*X2_3;

    R_or_C X3_12 = X3_1*X3_2 - X3_2*X3_1;
    R_or_C X3_21 = X3_2*X3_1 - X3_1*X3_2;
    R_or_C X3_13 = X3_1*X3_3 - X3_3*X3_1;
    R_or_C X3_31 = X3_3*X3_1 - X3_1*X3_3;
    R_or_C X3_23 = X3_2*X3_3 - X3_3*X3_2;
    R_or_C X3_32 = X3_3*X3_2 - X3_2*X3_3;

    matrix X4_12 = X4_1*X4_2 - X4_2*X4_1;
    matrix X4_21 = X4_2*X4_1 - X4_1*X4_2;
    matrix X4_13 = X4_1*X4_3 - X4_3*X4_1;
    matrix X4_31 = X4_3*X4_1 - X4_1*X4_3;
    matrix X4_23 = X4_2*X4_3 - X4_3*X4_2;
    matrix X4_32 = X4_3*X4_2 - X4_2*X4_3;

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
            A_vector_new[i] = Acceleration( i, X_vector_new, rows, cols, g);
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

        if (j % 1000 == 0)
        {
            //for (matrix el : V_vector)
            //{
            //    std::cout <<"Ideal " << el << std::endl;
            //}
            //std::cout  << std::endl << gauss_law(X_vector_new[0], X_vector_new[1], X_vector_new[2], X_vector_new[3], X_vector_new[4], X_vector_new[5], X_vector_new[6], X_vector_new[7], X_vector_new[8],
            //                       V_vector_new[0], V_vector_new[1], V_vector_new[2], V_vector_new[3], V_vector_new[4], V_vector_new[5], V_vector_new[6], V_vector_new[7], V_vector_new[8]);

            std::cout << std::endl;
            std::cout << "H" << std::setprecision(15) << H(g, 
                            X_vector_new[0], X_vector_new[1], X_vector_new[2], X_vector_new[3], X_vector_new[4], X_vector_new[5], X_vector_new[6], X_vector_new[7], X_vector_new[8],
                          V_vector_new[0], V_vector_new[1], V_vector_new[2], V_vector_new[3], V_vector_new[4], V_vector_new[5], V_vector_new[6], V_vector_new[7], V_vector_new[8]);


        }

    }




       // Export initial X/V/A_vector to text files to be analysed in python.
    std:: fstream X2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D2-Branes/thermalised_X.txt", std:: ios:: out);
    X2_vector_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (matrix Matrix : X_vector_new)
    {
        X2_vector_Export << Matrix << std::endl;
    }


    std:: fstream V2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D2-Branes/thermalised_V.txt", std:: ios:: out);
    V2_vector_Export << std::fixed << std::setprecision(15);
    for (matrix Matrix : V_vector_new)
    {
        V2_vector_Export << Matrix << std::endl;
    }

    std:: fstream A2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D2-Branes/thermalised_A.txt", std:: ios:: out);
    A2_vector_Export << std::fixed << std::setprecision(15);
    // Print to text file
    for (matrix Matrix : A_vector_new)
    {
        A2_vector_Export << Matrix << std::endl;
    }


    X2_vector_Export.close();
    V2_vector_Export.close();
    A2_vector_Export.close();


    std::cout << "Finished" << std::time(nullptr)-start;
    return 0;
}
