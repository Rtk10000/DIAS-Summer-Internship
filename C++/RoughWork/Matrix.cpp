#include <iostream>
#include <random>
#include <complex>

int main() {
    // Define the matrix dimensions
    const int rows = 3;
    const int cols = 3;
    const int size = rows * cols;

    // Create a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0.0, 1.0);

    // Declare and initialize the matrix
    std::complex<double> matrix[size];
    for (int i = 0; i < size; ++i) {
        double real = dist(gen);
        double imag = dist(gen);
        matrix[i] = std::complex<double>(real, imag);
    }

    // Print the matrix
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << matrix[i * cols + j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}