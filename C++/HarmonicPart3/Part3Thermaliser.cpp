/*To do list
    Define Hamilton's equations to get p_dot_all
    Find p_dot_initial
    Put p_dot_all_new code into iterator
    record every record_steps 
*/


#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

// Define timestep
const double delta_t = 10e-5;
// Repeat simulation for 1000 seconds.
const int simulation_repetitions = 1000 / delta_t;

// Define Lambda1 and Lambda2
const double lambda1 = 0.0001;
const double lambda2 = 0.0002;

// Define the Hamilton's equations
//Λ1𝑥32−3Λ1𝑥2𝑦22+3.0Λ22𝑥51+6.0Λ22𝑥31𝑦21+3.0Λ22𝑥1𝑦41+3Λ2𝑥21𝑥2+6Λ2𝑥1𝑦1𝑦2−3Λ2𝑥2𝑦21+𝑥1

double dHi_x_1(double x_1_old,double x_2_old,double y_1_old,double y_2_old,double lambda1,double lambda2)
{
    return lambda1 * pow(x_1_old, 3) - 3 * lambda1 * x_2_old * pow(y_2_old, 2)
    +3 * pow(lambda2, 2) * pow(x_1_old, 5) + 6 * pow(lambda2, 2) * pow(x_1_old, 3) * pow(y_1_old, 2)
    +3 * pow(lambda2, 2) * x_1_old * pow(y_2_old, 4) + 3 * lambda2 * pow(x_1_old, 2) * x_2_old
    +6 * pow(lambda2, 2) * x_1_old * y_1_old * y_2_old - 3 * lambda2 * x_2_old * pow(y_1_old, 2) + x_1_old;
}

int main()
{

    // Initialise and load in initial coordinates
    double x_1_old = 0, x_2_old = 0, y_1_old = 0, y_2_old = 0, p_x_1_old = 0, p_x_2_old = 0, p_y_1_old = 0, p_y_2_old = 0;
    
    std:: ifstream initial_coords("initial_conditions.txt", std::ios::in);
    initial_coords >> x_1_old >> x_2_old >> y_1_old >> y_2_old >> p_x_1_old >> p_x_2_old >> p_y_1_old >> p_y_2_old;
    initial_coords.close();

    // Evaluate initial rate of change of momentum
    double p_x_dot_1_old = -dHi_x_1(x_1_old, x_2_old, y_1_old, y_2_old, lambda1, lambda2);
    double p_x_dot_2_old = -dHi_x_2(x_1_old, x_2_old, y_1_old, y_2_old, lambda1, lambda2);
    double p_y_dot_1_old = -dHi_y_1(x_1_old, x_2_old, y_1_old, y_2_old, lambda1, lambda2);
    double p_y_dot_2_old = -dHi_y_2(x_1_old, x_2_old, y_1_old, y_2_old, lambda1, lambda2);

    // Initialise variables
    double x_1_new = 0, x_2_new = 0, y_1_new = 0, y_2_new = 0, p_x_1_new = 0, p_x_2_new = 0, p_y_1_new = 0, p_y_2_new = 0;

    // Run thermalising iterator
    for (int i = 1; i <= simulation_repetitions; i++) 
    {
        x_1_new = (x_1_old + p_x_1_old * delta_t + 1/2 * p_x_dot_1_old * delta_t**2)
        x_2_new = (x_2_old + p_x_2_old * delta_t + 1/2 * p_x_dot_2_old * delta_t**2)
        y_1_new = (y_1_old + p_y_1_old * delta_t + 1/2 * p_y_dot_1_old * delta_t**2)
        y_2_new = (y_2_old + p_y_2_old * delta_t + 1/2 * p_y_dot_2_old * delta_t**2)

        
        // c) Use Hamiltonian canonical equations to get new p_dot values
        
        p_x_dot_1_new = -dHi_x_1(x_1_new, x_2_new, y_1_new, y_2_new, Lambda_1, Lambda_2)
        p_x_dot_2_new = -dHi_x_2(x_1_new, x_2_new, y_1_new, y_2_new, Lambda_1, Lambda_2)
        p_y_dot_1_new = -dHi_y_1(x_1_new, x_2_new, y_1_new, y_2_new, Lambda_1, Lambda_2)
        p_y_dot_2_new = -dHi_y_2(x_1_new, x_2_new, y_1_new, y_2_new, Lambda_1, Lambda_2)


        // d) Use Velocity Verlet 2 to get new momentums

        p_x_1_new = (p_x_1_old + (1/2) * (p_x_dot_1_new + p_x_dot_1_old) * delta_t)
        p_x_2_new = (p_x_2_old + (1/2) * (p_x_dot_2_new + p_x_dot_2_old) * delta_t)
        p_y_1_new = (p_y_1_old + (1/2) * (p_y_dot_1_new + p_y_dot_1_old) * delta_t)
        p_y_2_new = (p_y_2_old + (1/2) * (p_y_dot_2_new + p_y_dot_2_old) * delta_t)

    }



    std:: fstream ThermalisationExport("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/ThermalisedPart3.txt", std:: ios:: out);
    ThermalisationExport << std::fixed << std::setprecision(15) << x_1_new << x_2_new << y_1_new << y_2_new << p_x_1_new << p_x_2_new << p_y_1_new << p_y_2_new << p_x_dot_1_new << p_x_dot_2_new << p_y_dot_1_new << p_y_dot_2_new;
    return 0;
}