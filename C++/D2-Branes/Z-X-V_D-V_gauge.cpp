//This is a starter to the full D2 branes Lagrangian. 
// Only using Zs as the generalised coords, and only having the potential term equal V_D


#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <random>
#include <complex>
#include <vector>

// The Eigen package, downloadable from the internet must be in the same folder as this .cpp file.
#include "eigen/Eigen/Dense"

// start variable is used to measure the amount of time it takes to run the code.
double start = std::time(nullptr);

// Random ci coefficients, I think these are valid numbers for them.
/*
long double c1 = gauss_dist(rng);
long double c2 = gauss_dist(rng);
long double c3 = gauss_dist(rng);
long double c4 = gauss_dist(rng);
*/

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

// Define timestep that each repetition will advance by.
const long double delta_t = 1e-4;

// Define the amount of seconds the system will be thermalised for
const long double seconds_thermalised = 10;

// Repeat simulation for appropriate amount of repetitions to reach the desired seconds thermalised.
const int simulation_repetitions = seconds_thermalised / delta_t;

// Number of D0-Branes
const int N = 3;

const long double g = 1.0;

// Define number of columns and rows each X matrix will have
const int rows = N;
const int cols = N;

// Dimension of space
const int dim = 9;

// Typedefs are used to define the type of object a variable is without writing the entire C++ name of each variable type

typedef std::complex<long double> Complex;
typedef Eigen:: Matrix<Complex, N, N> matrix;
typedef Eigen:: Matrix<Complex, 1, N> row_vector;
typedef Eigen:: Matrix<Complex, N, 1> col_vector;

// Initialise all Z generalised coords as randomised scalars or vectors
static std::random_device rd;
long double sigma = 1.0;
static std::mt19937 rng(10); 
std::normal_distribution<long double> dist(0, 1);

long double X12;
long double X21;

long double X13;
long double X31;

long double X23;
long double X32;

long double X12;
long double X21;




/*
Complex Z12 = Complex((dist(rng),dist(rng)));
Complex Z13 = Complex((dist(rng),dist(rng)));

Complex Z21 = Complex((dist(rng),dist(rng)));
Complex Z23 = Complex((dist(rng),dist(rng)));

Complex Z31 = Complex((dist(rng),dist(rng)));
Complex Z32 = Complex((dist(rng),dist(rng)));

// Random argument below gives the same Z vectors every time the code is run unfortunately.
row_vector Z14 = row_vector::Random();
row_vector Z24 = row_vector::Random();
row_vector Z34 = row_vector::Random();

col_vector Z41 = col_vector::Random();
col_vector Z42 = col_vector::Random();
col_vector Z43 = col_vector::Random();
*/

Complex Z12 = Complex(0.207606,0.366979);
Complex Z13 = Complex(0.0365052,-0.344526);

Complex Z21 = Complex(0.0786534,-1.71609);
Complex Z23 = Complex(0.508537,1.50702);

Complex Z31 = Complex(0.547891,0.496717);
Complex Z32 = Complex(-0.937337,-1.09565);
//Complex Z13 = Complex(0,0);
//Complex Z21 = Complex(0,0);
//Complex Z23 = Complex(0,0);
//Complex Z31 = Complex(0,0);
//Complex Z32 = Complex(0,0);

// Create a zero matrix and zero vectors.
matrix zero_matrix = matrix::Zero();
row_vector zero_row = row_vector::Zero();
col_vector zero_col = col_vector::Zero();

// Set the initial Z velocities to zeros
Complex Z12_dot = Complex(0.0, 0.0);
Complex Z13_dot = Complex(0.0, 0.0);
Complex Z21_dot = Complex(0.0, 0.0);
Complex Z23_dot = Complex(0.0, 0.0);
Complex Z31_dot = Complex(0.0, 0.0);
Complex Z32_dot = Complex(0.0, 0.0);

row_vector Z14_dot = zero_row;
row_vector Z24_dot = zero_row;
row_vector Z34_dot = zero_row;
col_vector Z41_dot = zero_col;
col_vector Z42_dot = zero_col;
col_vector Z43_dot = zero_col;



// Define H = K + V
long double H(
    const long double g, 
    Complex Z12, Complex Z13, Complex Z21, Complex Z23, Complex Z31, Complex Z32,
    row_vector Z14, row_vector Z24, row_vector Z34,
    col_vector Z41, col_vector Z42, col_vector Z43,
    Complex Z12_dot, Complex Z13_dot, Complex Z21_dot, Complex Z23_dot, Complex Z31_dot, Complex Z32_dot,
    row_vector Z14_dot, row_vector Z24_dot, row_vector Z34_dot, 
    col_vector Z41_dot, col_vector Z42_dot, col_vector Z43_dot)
{   

    long double K = 0.5 * (std::conj(Z12_dot)*Z12_dot + std::conj(Z13_dot)*Z13_dot +
                                   std::conj(Z21_dot)*Z21_dot + std::conj(Z23_dot)*Z23_dot + 
                                   std::conj(Z31_dot)*Z31_dot + std::conj(Z32_dot)*Z32_dot + 
                                   (Z14_dot.adjoint()*Z14_dot).trace().real() + (Z24_dot.adjoint()*Z24_dot).trace().real() + 
                                   (Z34_dot.adjoint()*Z34_dot).trace().real() + (Z41_dot.adjoint()*Z41_dot).trace().real() + 
                                   (Z42_dot.adjoint()*Z42_dot).trace().real() + (Z43_dot.adjoint()*Z43_dot).trace().real() ).real();

    // for k = 1
    Complex V_D_arg_k1 = (Z12 * std::conj(Z12) + Z13 * std::conj(Z13) + (Z14 * Z14.adjoint()).trace() -
                         (std::conj(Z21) * Z21 + std::conj(Z31) * Z31 + (Z41.adjoint() * Z41).trace())
                         - c1/(g*g));
    // for k = 2
    Complex V_D_arg_k2 = (Z21 * std::conj(Z21) + Z23 * std::conj(Z23) + (Z24 * Z24.adjoint()).trace() - 
                         (std::conj(Z12) * Z12 + std::conj(Z32) * Z32 + (Z42.adjoint() * Z42).trace() )
                         - c2/(g*g) );
    // for k = 3
    Complex V_D_arg_k3 = (Z31 * std::conj(Z31) + Z32 * std::conj(Z32) + (Z34 * Z34.adjoint()).trace() - 
                         (std::conj(Z13) * Z13 + std::conj(Z23) * Z23 + (Z43.adjoint() * Z43).trace() ) 
                         - c3/(g*g) );
    // for k = 4
    matrix V_D_arg_k4 = (Z41 * Z41.adjoint() + Z42 * Z42.adjoint() + Z43 * Z43.adjoint() - 
                        (Z14.adjoint() * Z14 + Z24.adjoint() * Z24 + Z34.adjoint() * Z34 ) 
                        - c4/(g*g) * matrix::Identity() );

    long double V_D = 0.5 * ((V_D_arg_k1 * V_D_arg_k1) + (V_D_arg_k2 * V_D_arg_k2) + 
             (V_D_arg_k3 * V_D_arg_k3) + (V_D_arg_k4 * V_D_arg_k4).trace().real() ).real();

    long double V = g*g*V_D;

    return K + V;
}

int main() 
{

    row_vector Z14;
    Z14(0,0) = Complex(1,1);
    Z14(0,1) = Complex(2,1);
    Z14(0,2) = Complex(3,1);
    row_vector Z24;
    Z24(0,0) = Complex(1,2);
    Z24(0,1) = Complex(2,2);
    Z24(0,2) = Complex(3,2);
    row_vector Z34;
    Z34(0,0) = Complex(1,3);
    Z34(0,1) = Complex(2,3);
    Z34(0,2) = Complex(3,3);

    col_vector Z41;
    Z41(0,0) = Complex(1,1);
    Z41(1,0) = Complex(1,2);
    Z41(2,0) = Complex(1,3);
    col_vector Z42;
    Z42(0,0) = Complex(2,1);
    Z42(1,0) = Complex(2,2);
    Z42(2,0) = Complex(2,3);
    col_vector Z43;
    Z43(0,0) = Complex(3,1);
    Z43(1,0) = Complex(3,2);
    Z43(2,0) = Complex(3,3);

/*
    row_vector Z14;
    Z14(0,0) = Complex(0,0);
    Z14(0,1) = Complex(0,0);
    Z14(0,2) = Complex(0,0);
    row_vector Z24;
    Z24(0,0) = Complex(0,0);
    Z24(0,1) = Complex(0,0);
    Z24(0,2) = Complex(0,0);
    row_vector Z34;
    Z34(0,0) = Complex(0,0);
    Z34(0,1) = Complex(0,0);
    Z34(0,2) = Complex(0,0);

    col_vector Z41;
    Z41(0,0) = Complex(0,0);
    Z41(1,0) = Complex(0,0);
    Z41(2,0) = Complex(0,0);
    col_vector Z42;
    Z42(0,0) = Complex(0,0);
    Z42(1,0) = Complex(0,0);
    Z42(2,0) = Complex(0,0);
    col_vector Z43;
    Z43(0,0) = Complex(0,0);
    Z43(1,0) = Complex(0,0);
    Z43(2,0) = Complex(0,0);
*/
    // These terms are present in all the equations of motion for Z.
    Complex c1_term = abs(Z12)*abs(Z12) + abs(Z13)*abs(Z13) + (Z14 * Z14.adjoint()).trace() - abs(Z21)*abs(Z21) - abs(Z31)*abs(Z31) - (Z41.adjoint() * Z41).trace() - c1/(g*g);
    Complex c2_term = abs(Z21)*abs(Z21) + abs(Z23)*abs(Z23) + (Z24 * Z24.adjoint()).trace() - abs(Z12)*abs(Z12) - abs(Z32)*abs(Z32) - (Z42.adjoint() * Z42).trace() - c2/(g*g);
    Complex c3_term = abs(Z31)*abs(Z31) + abs(Z32)*abs(Z32) + (Z34 * Z34.adjoint()).trace() - abs(Z13)*abs(Z13) - abs(Z23)*abs(Z23) - (Z43.adjoint() * Z43).trace() - c3/(g*g);
    matrix  c4_term = (Z41 * Z41.adjoint()) + (Z42 * Z42.adjoint()) + (Z43 * Z43.adjoint()) - (Z14.adjoint() * Z14) - (Z24.adjoint() * Z24) - (Z34.adjoint() * Z34) - c4/(g*g) * matrix::Identity();

    // Define the Z accelerations from the equations of motion.
    Complex Z12_double_dot = Complex(2,0) * (-Z12 * c1_term + Z12 * c2_term);
    Complex Z13_double_dot = Complex(2,0) * (-Z13 * c1_term + Z13 * c3_term);
    Complex Z21_double_dot = Complex(2,0) * (+Z21 * c1_term - Z21 * c2_term);
    Complex Z23_double_dot = Complex(2,0) * (-Z23 * c2_term + Z23 * c3_term);
    Complex Z31_double_dot = Complex(2,0) * (+Z31 * c1_term - Z31 * c3_term);
    Complex Z32_double_dot = Complex(2,0) * (+Z32 * c2_term - Z32 * c3_term);

    row_vector Z14_double_dot = Complex(2,0) * (-Z14 * c1_term + Z14 * c4_term);
    row_vector Z41_double_dot = Complex(2,0) * ( Z41 * c1_term - c4_term * Z41);
    row_vector Z24_double_dot = Complex(2,0) * (-Z24 * c2_term + Z24 * c4_term);
    col_vector Z42_double_dot = Complex(2,0) * ( Z42 * c1_term - c4_term * Z42);
    col_vector Z34_double_dot = Complex(2,0) * (-Z34 * c3_term + Z34 * c4_term);
    col_vector Z43_double_dot = Complex(2,0) * ( Z43 * c3_term - c4_term * Z43);

 // Put all variables into arrays to more easily write loops for the numerical integration algorithm.

    Complex scalar_Zs[6] = {Z12, Z13, Z21, Z23, Z31, Z32};
    row_vector row_Zs[3] = {Z14, Z24, Z34};
    col_vector col_Zs[3] = {Z41, Z42, Z43};

 // Initialise the new Z arrays which will then be populated by the new Z values later.
    Complex scalar_Zs_new[6] = {Complex(0.0, 0.0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)};
    row_vector row_Zs_new[3] = {zero_row, zero_row, zero_row};
    col_vector col_Zs_new[3] = {zero_col, zero_col, zero_col};

    Complex scalar_Z_dots[6] = {Z12_dot, Z13_dot, Z21_dot, Z23_dot, Z31_dot, Z32_dot};
    row_vector row_Z_dots[3] = {Z14_dot, Z24_dot, Z34_dot};
    col_vector col_Z_dots[3] = {Z41_dot, Z42_dot, Z43_dot};

    Complex scalar_Z_dots_new[6] = {Complex(0.0, 0.0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)};
    row_vector row_Z_dots_new[3] = {zero_row, zero_row, zero_row};
    col_vector col_Z_dots_new[3] = {zero_col, zero_col, zero_col};

    Complex scalar_Z_double_dots[6] = {Z12_double_dot, Z13_double_dot, Z21_double_dot, Z23_double_dot, Z31_double_dot, Z32_double_dot};
    row_vector row_Z_double_dots[3] = {Z14_double_dot, Z24_double_dot, Z34_double_dot};
    col_vector col_Z_double_dots[3] = {Z41_double_dot, Z42_double_dot, Z43_double_dot};

    Complex scalar_Z_double_dots_new[6] = {Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0)};
    row_vector row_Z_double_dots_new[3] = {zero_row, zero_row, zero_row};
    col_vector col_Z_double_dots_new[3] = {zero_col, zero_col, zero_col};



// Create links to text files which we will populate with the Z and Z_dot and Z_double_dot values during the simulation.
    /*    
        // Export initial Z, Z_dot, and Z_double_dot objects to text files to be analysed in python.
    std:: fstream Z_scalar_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/D2-Branes/Z_scalar.txt", std:: ios:: out);
    Z_scalar_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (Complex Z : scalar_Zs)
    {
        Z_scalar_Export << Z << std::endl;
    }

    std:: fstream Z_row_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/D2-Branes/Z_row.txt", std:: ios:: out);
    Z_row_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (row_vector Z : row_Zs)
    {
        Z_row_Export << Z << std::endl;
    }

    std:: fstream Z_col_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/D2-Branes/Z_col.txt", std:: ios:: out);
    Z_col_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (col_vector Z : col_Zs)
    {
        Z_col_Export << Z << std::endl;
    }
    
    std:: fstream Z_dot_scalar_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/D2-Branes/Z_dot_scalar.txt", std:: ios:: out);
    Z_dot_scalar_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (Complex Z : scalar_Z_dots)
    {
        Z_dot_scalar_Export << Z << std::endl;
    }

    std:: fstream Z_dot_row_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/D2-Branes/Z_dot_row.txt", std:: ios:: out);
    Z_dot_row_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (row_vector Z : row_Z_dots)
    {
        Z_dot_row_Export << Z << std::endl;
    }

    std:: fstream Z_dot_col_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/D2-Branes/Z_dot_col.txt", std:: ios:: out);
    Z_dot_col_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (col_vector Z : col_Z_dots)
    {
        Z_dot_col_Export << Z << std::endl;
    }


    std:: fstream Z_double_dot_scalar_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/D2-Branes/Z_double_dot_scalar.txt", std:: ios:: out);
    Z_double_dot_scalar_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (Complex Z : scalar_Z_double_dots)
    {
        Z_double_dot_scalar_Export << Z << std::endl;
    }

    std:: fstream Z_double_dot_row_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/D2-Branes/Z_double_dot_row.txt", std:: ios:: out);
    Z_double_dot_row_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (row_vector Z : row_Z_double_dots)
    {
        Z_double_dot_row_Export << Z << std::endl;
    }

    std:: fstream Z_double_dot_col_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/D2-Branes/Z_double_dot_col.txt", std:: ios:: out);
    Z_double_dot_col_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (col_vector Z : col_Z_double_dots)
    {
        Z_double_dot_col_Export << Z << std::endl;
    }
    */

    std::cout << std::setprecision(15)<< "Starting Energy" << H(g, 
                            scalar_Zs[0], scalar_Zs[1], scalar_Zs[2], scalar_Zs[3], scalar_Zs[4], scalar_Zs[5],
                            row_Zs[0], row_Zs[1], row_Zs[2],
                            col_Zs[0], col_Zs[1], col_Zs[2],
                            scalar_Z_dots[0], scalar_Z_dots[1], scalar_Z_dots[2], scalar_Z_dots[3], scalar_Z_dots[4], scalar_Z_dots[5],
                            row_Z_dots[0], row_Z_dots[1], row_Z_dots[2],
                            col_Z_dots[0], col_Z_dots[1], col_Z_dots[2]);
    std::cout << '\n';


    // Write simulation to thermalise system
    for (int j = 0; j < seconds_thermalised / delta_t; ++j)
    {


            // Every 1000 steps calculate the total energy of the system to see if it is conserved.
        if (j % 1000 == 0)
        {
            std::cout << '\n' << "H at" << j*delta_t << "seconds: ";
            std::cout << std::setprecision(15) << H(g, 
                            scalar_Zs[0], scalar_Zs[1], scalar_Zs[2], scalar_Zs[3], scalar_Zs[4], scalar_Zs[5],
                            row_Zs[0], row_Zs[1], row_Zs[2],
                            col_Zs[0], col_Zs[1], col_Zs[2],
                            scalar_Z_dots[0], scalar_Z_dots[1], scalar_Z_dots[2], scalar_Z_dots[3], scalar_Z_dots[4], scalar_Z_dots[5],
                            row_Z_dots[0], row_Z_dots[1], row_Z_dots[2],
                            col_Z_dots[0], col_Z_dots[1], col_Z_dots[2]) << '\n';
            std::cout << "Z12 = " << scalar_Zs[0] <<'\n';

        }
        /*
        Algorithm works as follows:
        Step 1) Calculate Z_new from old Z, Z_dot, and Z_double_dot values from the first velocity Verlet equation

        Step 2) Calculate Z_double_dot_new from the Z_new calculated in step 1)

        Step 3) Calculate Z_dot_new from Z_dot, Z_double_dot, and Z_double_dot_new.

        Step 4) Redefine Z_new as Z, Z_dot_new as Z_dot, and Z_double_dot_new as Z_double_dot.

        Step 5) Repeat all steps
        */

        // Step 1) velocity Verlet 1 to get new positions from old positions, momentums and rate of change of momentums

        for (int i = 0; i < 6; ++i)
        {
            scalar_Zs_new[i] = scalar_Zs[i] + scalar_Z_dots[i] * delta_t + Complex(0.5, 0) * scalar_Z_double_dots[i] * delta_t * delta_t;
        }
        for (int i = 0; i < 3; ++i)
        {
            row_Zs_new[i] = row_Zs[i] + row_Z_dots[i] * delta_t + Complex(0.5, 0) * row_Z_double_dots[i] * delta_t * delta_t;
        }
        for (int i = 0; i < 3; ++i)
        {
            col_Zs_new[i] = col_Zs[i] + col_Z_dots[i] * delta_t + Complex(0.5, 0) * col_Z_double_dots[i] * delta_t * delta_t;
        }

        // Step 2) Generate and store new Z_double_dots

        Complex c1_term_new = abs(scalar_Zs_new[0])*abs(scalar_Zs_new[0]) + abs(scalar_Zs_new[1])*abs(scalar_Zs_new[1]) + (row_Zs_new[0] * row_Zs_new[0].adjoint()).trace() - abs(scalar_Zs_new[2])*abs(scalar_Zs_new[2]) - abs(scalar_Zs_new[4])*abs(scalar_Zs_new[4]) - (col_Zs_new[0].adjoint() * col_Zs_new[0]).trace() - c1/(g*g);
        Complex c2_term_new = abs(scalar_Zs_new[2])*abs(scalar_Zs_new[2]) + abs(scalar_Zs_new[3])*abs(scalar_Zs_new[3]) + (row_Zs_new[1] * row_Zs_new[1].adjoint()).trace() - abs(scalar_Zs_new[0])*abs(scalar_Zs_new[0]) - abs(scalar_Zs_new[5])*abs(scalar_Zs_new[5]) - (col_Zs_new[1].adjoint() * col_Zs_new[1]).trace() - c2/(g*g);
        Complex c3_term_new = abs(scalar_Zs_new[4])*abs(scalar_Zs_new[4]) + abs(scalar_Zs_new[5])*abs(scalar_Zs_new[5]) + (row_Zs_new[2] * row_Zs_new[2].adjoint()).trace() - abs(scalar_Zs_new[1])*abs(scalar_Zs_new[1]) - abs(scalar_Zs_new[3])*abs(scalar_Zs_new[3]) - (col_Zs_new[2].adjoint() * col_Zs_new[2]).trace() - c3/(g*g);
        matrix  c4_term_new = (col_Zs_new[0] * col_Zs_new[0].adjoint()) + (col_Zs_new[1] * col_Zs_new[1].adjoint()) + (col_Zs_new[2] * col_Zs_new[2].adjoint()) - (row_Zs_new[0].adjoint() * row_Zs_new[0]) - (row_Zs_new[1].adjoint() * row_Zs_new[1]) - (row_Zs_new[2].adjoint() * row_Zs_new[2]) - c4/(g*g) * matrix::Identity();

        scalar_Z_double_dots_new[0] = Complex(2,0) * (-scalar_Zs_new[0] * c1_term_new + scalar_Zs_new[0] * c2_term_new);
        scalar_Z_double_dots_new[1] = Complex(2,0) * (-scalar_Zs_new[1] * c1_term_new + scalar_Zs_new[1] * c3_term_new);
        scalar_Z_double_dots_new[2] = Complex(2,0) * (+scalar_Zs_new[2] * c1_term_new - scalar_Zs_new[2] * c2_term_new);
        scalar_Z_double_dots_new[3] = Complex(2,0) * (-scalar_Zs_new[3] * c2_term_new + scalar_Zs_new[3] * c3_term_new);
        scalar_Z_double_dots_new[4] = Complex(2,0) * (+scalar_Zs_new[4] * c1_term_new - scalar_Zs_new[4] * c3_term_new);
        scalar_Z_double_dots_new[5] = Complex(2,0) * (+scalar_Zs_new[5] * c2_term_new - scalar_Zs_new[5] * c3_term_new);

        row_Z_double_dots_new[0] = Complex(2,0) * (-row_Zs_new[0] * c1_term_new + row_Zs_new[0] * c4_term_new);
        col_Z_double_dots_new[0] = Complex(2,0) * ( col_Zs_new[0] * c1_term_new - c4_term_new * col_Zs_new[0]);

        row_Z_double_dots_new[1] = Complex(2,0) * (-row_Zs_new[1] * c2_term_new + row_Zs_new[1] * c4_term_new);
        col_Z_double_dots_new[1] = Complex(2,0) * ( col_Zs_new[1] * c2_term_new - c4_term_new * col_Zs_new[1]);

        row_Z_double_dots_new[2] = Complex(2,0) * (-row_Zs_new[2] * c3_term_new + row_Zs_new[2] * c4_term_new);
        col_Z_double_dots_new[2] = Complex(2,0) * ( col_Zs_new[2] * c3_term_new - c4_term_new * col_Zs_new[2]);  


        // Step 3) Use Velocity Verlet 2 to get new momentums
        for (int i = 0; i < 6; ++i)
        {
            scalar_Z_dots_new[i] = scalar_Z_dots[i] + Complex(0.5, 0) * (scalar_Z_double_dots_new[i] + scalar_Z_double_dots[i]) * delta_t;
        }
        for (int i = 0; i < 3; ++i)
        {
            row_Z_dots_new[i] = row_Z_dots[i] + 0.5 * (row_Z_double_dots_new[i] + row_Z_double_dots[i]) * delta_t;
        }
        for (int i = 0; i < 3; ++i)
        {
            col_Z_dots_new[i] = col_Z_dots[i] + 0.5 * (col_Z_double_dots_new[i] + col_Z_double_dots[i]) * delta_t;
        }


        // Step 4)
        // Copy elements from X_vector_new to X_vector
        std::memcpy(scalar_Zs, scalar_Zs_new, sizeof(scalar_Zs_new));  
        std::memcpy(row_Zs, row_Zs_new, sizeof(row_Zs_new));  
        std::memcpy(col_Zs, col_Zs_new, sizeof(col_Zs_new));  

        // Copy elements from V_vector_new to V_vector
        std::memcpy(scalar_Z_dots, scalar_Z_dots_new, sizeof(scalar_Z_dots_new));  
        std::memcpy(row_Z_dots, row_Z_dots_new, sizeof(row_Z_dots_new));  
        std::memcpy(col_Z_dots, col_Z_dots_new, sizeof(col_Z_dots_new));  

        // Copy elements from A_vector_new to A_vector
        std::memcpy(scalar_Z_double_dots, scalar_Z_double_dots_new, sizeof(scalar_Z_double_dots_new));  
        std::memcpy(row_Z_double_dots, row_Z_double_dots_new, sizeof(row_Z_double_dots_new));  
        std::memcpy(col_Z_double_dots, col_Z_double_dots_new, sizeof(col_Z_double_dots_new));  

        // Every 1000 steps Print Zs to text files for analysis in Python.
        /*
        if (j % 1000 == 0)
        {
        // Print Zs to text files
            //Print to text file
            for (Complex Z : scalar_Zs)
            {
            Z_scalar_Export << Z << std::endl;
            }

            //Print to text file
            for (row_vector Z : row_Zs)
            {
            Z_row_Export << Z << std::endl;
            }

            //Print to text file
            for (col_vector Z : col_Zs)
            {
            Z_col_Export << Z << std::endl;
            }
    
            //Print to text file
            for (Complex Z : scalar_Z_dots)
            {
                Z_dot_scalar_Export << Z << std::endl;
            }

            //Print to text file
            for (row_vector Z : row_Z_dots)
            {
                Z_dot_row_Export << Z << std::endl;
            }

            //Print to text file
            for (col_vector Z : col_Z_dots)
            {
                Z_dot_col_Export << Z << std::endl;
            }


                //Print to text file
            for (Complex Z : scalar_Z_double_dots)
            {
                Z_double_dot_scalar_Export << Z << std::endl;
            }

            //Print to text file
            for (row_vector Z : row_Z_double_dots)
            {
                Z_double_dot_row_Export << Z << std::endl;
            }

            //Print to text file
            for (col_vector Z : col_Z_double_dots)
            {
                Z_double_dot_col_Export << Z << std::endl;
            }
        }
        */
        // Step 5) repeat all steps.
    }


    // Zs have now been put into text files, so they can be uploaded to python or to other.cpp files.
    /*
    Z_scalar_Export.close();
    Z_row_Export.close();
    Z_col_Export.close();
    Z_dot_scalar_Export.close();
    Z_dot_row_Export.close();
    Z_dot_col_Export.close();
    Z_double_dot_scalar_Export.close();
    Z_double_dot_row_Export.close();
    Z_double_dot_col_Export.close();
    */
    // Print the time it took to run this .cpp file.
    std::cout << "Finished" << std::time(nullptr)-start;

    return 0;
}