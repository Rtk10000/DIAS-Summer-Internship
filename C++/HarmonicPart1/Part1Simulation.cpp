#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

// Define timestep
const double delta_t = 1e-3;
const int simulation_repetitions = 1000000;
const int record_steps = 1000;

int main() 
{

    // Initialise and load in initial coordinates
    double x_old = 0, p_old = 0;
    const int k = 1;
    
    std:: ifstream initial_coords("initial_coords.txt", std::ios::in);
    initial_coords >> x_old >> p_old;
    initial_coords.close();
    
    // Evaluate initial rate of change of momentum
    double p_dot_old = -k * x_old;

    // Initialise variables
    double x_new = 0, p_new = 0, p_dot_new = 0;

    // Open files so x and p can be put in during the for loop simulation.
    std:: fstream xSimulationExport("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/x1WarmupSimulation.txt", std:: ios:: out);
    xSimulationExport << std::fixed << std::setprecision(15) << x_old;

    std:: fstream pSimulationExport("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/p1WarmupSimulation.txt", std:: ios:: out);
    pSimulationExport << std::fixed << std::setprecision(15) << p_old;

    // Evolving the system in time

    for (int i = 1; i <= simulation_repetitions; i++) 
    {
        x_new = x_old + p_old * delta_t + 0.5 * p_dot_old * pow(delta_t, 2); // Velocity Verlet algorithm 1.

        p_dot_new = - k * x_new;

        p_new = p_old + 0.5 * (p_dot_new + p_dot_old) * delta_t; // Velocity Verlet alogorithm 2.

    
        // Record every 1000th set of coordinates
        if (i % record_steps == 0)
        {
            // Export new x and p into text files.
            xSimulationExport << "," << x_new;
            pSimulationExport << "," << p_new;
            
        }

        // Update coords
        x_old = x_new;
        p_old = p_new;
        p_dot_old = p_dot_new;
    }

    
    
    
    
    // Close file streams.
    xSimulationExport.close();
    pSimulationExport.close();

    return 0;
}