#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <random>
#include <complex>
#include <vector>

const double delta_t = 1e-6;

int main() {
    std:: cout << -10 * pow(delta_t, 2);

    std:: cout << "Enter a radius:";

    const double pi = 3.14159;
    double radius = 0;

    //std:: cin >> radius;

    std:: cout << "area of circle is: " << pi * pow(radius, 2);




    double f[2] = {2,4};
    double g[2] = {4,7};
    std::cout << "f:" << f<< std::endl;
    std::cout << "g before:" << g<< std::endl;

    std::memcpy(&g, &f, sizeof(f));  
    std::cout << " f after:" << g<< std::endl;
    
    std::cout << "g after:" << g<< std::endl;



    return 0;
}