#include <iostream>
#include "eigen/Eigen/Dense"


const int N = 3;


void printArray(int* arr, int size) 
{   
    for (int i = 0; i < size; i++) 
    {
        int X = arr[i];
        std::cout << X << " ";
    }
    std::cout << std::endl;
}
typedef double R_or_C;
typedef Eigen:: Matrix<double, N, N> matrix;


int main() {

    const int rows = N;
    const int cols = N;

    int myArray[] = {1, 2, 99, 4, 5};
    int size = sizeof(myArray) / sizeof(myArray[0]);

    printArray(myArray, size);



    matrix zero;

    std::cout << zero;

    return 0;
}