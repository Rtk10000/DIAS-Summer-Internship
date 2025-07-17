#include <iostream>
#include <fstream>

int main() {

    // Load in initial coordinates
    std:: ifstream initial_coords("initial_coords.txt", std::ios::in);

    // Initialise the position and momentum variables
    double x1, p1;

    initial_coords >> x1 >> p1; // Read values from file
    initial_coords.close();     // Close file after reading in values.

    // Print out other coordinates into txt file.
    std::cout << x1 << p1;

    double x1_new = x1+2;
    double p1_new = p1+2;


    // Export new x1 into text file.
    std:: fstream xSimulationExport("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/x1WarmupSimulation.txt", std:: ios:: out);
    xSimulationExport << x1_new;
    xSimulationExport.close();

    // Export new p1 into text file.
    std:: fstream pSimulationExport("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/p1WarmupSimulation.txt", std:: ios:: out);
    pSimulationExport << p1_new;
    pSimulationExport.close();


    return 0;
}