#include <iostream>
#include <random>
#include <ctime>
#include <cmath>
#include <fstream>
#include "eigen/Eigen/Dense"
#include <iomanip>

int start = std::time(nullptr);
const int N = 4;
const int iterations = 2e6;
const long double dt = 1e-4;
const int perturb_iterations = 1e4; 
const int record = 1000;

typedef std::complex<long double> complex;
typedef Eigen:: Matrix<std::complex<long double>, N, N> matrix;

std:: mt19937 rng(std:: time(nullptr));
std:: normal_distribution<long double> coeff_dist(0, 1e-8);
const long double c1 = 2e-9;
const long double c2 = 5e-9;
long double coefficients[2] = {c1,c2};

matrix commutator(matrix A, matrix B)
{
    return A*B - B*A;
} 

matrix anticommutator(matrix A, matrix B)
{
    return A*B + B*A;
}

long double H(
    long double g, 
    matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7, matrix X8, matrix X9,
    matrix V1, matrix V2, matrix V3, matrix V4, matrix V5, matrix V6, matrix V7, matrix V8, matrix V9)
{
    // Compute kinetic energy T
    complex T = 1/(2*g*g) * (V1*V1 + V2*V2 + V3*V3 + V4*V4 + V5*V5 + V6*V6 + V7*V7 + V8*V8 + V9*V9).trace();

    matrix X[9] = {X1,X2,X3,X4,X5,X6,X7,X8,X9}; 

    matrix commutator_sum;  
    for (int i = 0; i < 9; i++)
    {
        for (int j = 0; j < 9; j++)
        {
            if(i == j)
                continue;
            commutator_sum += commutator(X[i],X[j])*commutator(X[i],X[j]); //can likely be more efficient by less function calls
        }
    }
    complex U = - (1)/(4*g*g) * commutator_sum.trace();
    return std:: abs(T + U);
}

matrix F(int i, matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7,matrix X8, matrix X9)
{
    matrix X[9] = {X1,X2,X3,X4,X5,X6,X7,X8,X9};
    matrix sum;
    for (int k = 0; k < 9; k++)
    {
        if(i-1 == k)
        {
            continue;
        }
        sum += commutator(X[k],commutator(X[i - 1],X[k]));
    }
    return sum;
}

matrix DF(int i, matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7, matrix X8, matrix X9)
{
    matrix X[9] = {X1,X2,X3,X4,X5,X6,X7,X8,X9};
    matrix sum1, sum2;

    for (int k = 0; k < 9; k++)
    {
        if(i-1 == k)
        {
            continue;
        }
        sum1 += commutator(X[k],commutator(X[i - 1],X[k]));
    }
    
    matrix product_sum;
    for (int j = 0; j < 9; j++)
    {
        product_sum += X[j]*X[j];
    }

    sum2 = c1*2*X[i-1] + 2*c2*anticommutator(X[i-1],product_sum);

    return (sum1 + sum2);
}

matrix gauss_law(
    matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7, matrix X8, matrix X9,
    matrix V1, matrix V2, matrix V3, matrix V4, matrix V5, matrix V6, matrix V7, matrix V8, matrix V9)
{
    matrix result;
    for (int i = 0; i < 9; i ++)
    {
        result = commutator(X1,V1) + commutator(X2,V2) + commutator(X3,V3) + commutator(X4,V4) + commutator(X5,V5) +
        commutator(X6,V6) + commutator(X7,V7) + commutator(X8,V8) + commutator(X9,V9);
    }
    return result;
}

void update(
    long double dt,
    matrix* X1, matrix* X2, matrix* X3, matrix* X4, matrix* X5, matrix* X6, matrix* X7, matrix* X8, matrix* X9,
    matrix* X1_n, matrix* X2_n, matrix* X3_n, matrix* X4_n, matrix* X5_n, matrix* X6_n, matrix* X7_n, matrix* X8_n, matrix* X9_n,
    matrix* V1, matrix* V2, matrix* V3, matrix* V4, matrix* V5, matrix* V6, matrix* V7, matrix* V8, matrix* V9,
    matrix* F1_0, matrix* F2_0, matrix* F3_0, matrix* F4_0, matrix* F5_0, matrix* F6_0, matrix* F7_0, matrix* F8_0, matrix* F9_0,
    matrix* F1_n, matrix* F2_n, matrix* F3_n, matrix* F4_n, matrix* F5_n, matrix* F6_n, matrix* F7_n, matrix* F8_n, matrix* F9_n)
{
    *X1_n = *X1 + (*V1)*dt + 0.5*(*F1_0)*pow(dt,2);
    *X2_n = *X2 + (*V2)*dt + 0.5*(*F2_0)*pow(dt,2);
    *X3_n = *X3 + (*V3)*dt + 0.5*(*F3_0)*pow(dt,2);
    *X4_n = *X4 + (*V4)*dt + 0.5*(*F4_0)*pow(dt,2);
    *X5_n = *X5 + (*V5)*dt + 0.5*(*F5_0)*pow(dt,2);
    *X6_n = *X6 + (*V6)*dt + 0.5*(*F6_0)*pow(dt,2);
    *X7_n = *X7 + (*V7)*dt + 0.5*(*F7_0)*pow(dt,2);
    *X8_n = *X8 + (*V8)*dt + 0.5*(*F8_0)*pow(dt,2);
    *X9_n = *X9 + (*V9)*dt + 0.5*(*F9_0)*pow(dt,2);

    *F1_n = F(1,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F2_n = F(2,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F3_n = F(3,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F4_n = F(4,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F5_n = F(5,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F6_n = F(6,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F7_n = F(7,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F8_n = F(8,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F9_n = F(9,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);

    *V1 = *V1 + 0.5*(*F1_0 + *F1_n)*dt;
    *V2 = *V2 + 0.5*(*F2_0 + *F2_n)*dt;
    *V3 = *V3 + 0.5*(*F3_0 + *F3_n)*dt;
    *V4 = *V4 + 0.5*(*F4_0 + *F4_n)*dt;
    *V5 = *V5 + 0.5*(*F5_0 + *F5_n)*dt;
    *V6 = *V6 + 0.5*(*F6_0 + *F6_n)*dt;
    *V7 = *V7 + 0.5*(*F7_0 + *F7_n)*dt;
    *V8 = *V8 + 0.5*(*F8_0 + *F8_n)*dt;
    *V9 = *V9 + 0.5*(*F9_0 + *F9_n)*dt;

    *X1 = *X1_n, *X2 = *X2_n, *X3 = *X3_n, *X4 = *X4_n, *X5 = *X5_n, *X6 = *X6_n, *X7 = *X7_n, *X8 = *X8_n, *X9 = *X9_n;
    *F1_0 = *F1_n, *F2_0 = *F2_n, *F3_0 = *F3_n, *F4_0 = *F4_n, *F5_0 = *F5_n, *F6_0 = *F6_n, *F7_0 = *F7_n, *F8_0 = *F8_n, *F9_0 = *F9_n;
}

void perturb_update(
    long double dt,
    matrix* X1, matrix* X2, matrix* X3, matrix* X4, matrix* X5, matrix* X6, matrix* X7, matrix* X8, matrix* X9,
    matrix* X1_n, matrix* X2_n, matrix* X3_n, matrix* X4_n, matrix* X5_n, matrix* X6_n, matrix* X7_n, matrix* X8_n, matrix* X9_n,
    matrix* V1, matrix* V2, matrix* V3, matrix* V4, matrix* V5, matrix* V6, matrix* V7, matrix* V8, matrix* V9,
    matrix* F1_0, matrix* F2_0, matrix* F3_0, matrix* F4_0, matrix* F5_0, matrix* F6_0, matrix* F7_0, matrix* F8_0, matrix* F9_0,
    matrix* F1_n, matrix* F2_n, matrix* F3_n, matrix* F4_n, matrix* F5_n, matrix* F6_n, matrix* F7_n, matrix* F8_n, matrix* F9_n)
{
    *X1_n = *X1 + (*V1)*dt + 0.5*(*F1_0)*pow(dt,2);
    *X2_n = *X2 + (*V2)*dt + 0.5*(*F2_0)*pow(dt,2);
    *X3_n = *X3 + (*V3)*dt + 0.5*(*F3_0)*pow(dt,2);
    *X4_n = *X4 + (*V4)*dt + 0.5*(*F4_0)*pow(dt,2);
    *X5_n = *X5 + (*V5)*dt + 0.5*(*F5_0)*pow(dt,2);
    *X6_n = *X6 + (*V6)*dt + 0.5*(*F6_0)*pow(dt,2);
    *X7_n = *X7 + (*V7)*dt + 0.5*(*F7_0)*pow(dt,2);
    *X8_n = *X8 + (*V8)*dt + 0.5*(*F8_0)*pow(dt,2);
    *X9_n = *X9 + (*V9)*dt + 0.5*(*F9_0)*pow(dt,2);

    *F1_n = DF(1,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F2_n = DF(2,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F3_n = DF(3,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F4_n = DF(4,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F5_n = DF(5,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F6_n = DF(6,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F7_n = DF(7,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F8_n = DF(8,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F9_n = DF(9,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);

    *V1 = *V1 + 0.5*(*F1_0 + *F1_n)*dt;
    *V2 = *V2 + 0.5*(*F2_0 + *F2_n)*dt;
    *V3 = *V3 + 0.5*(*F3_0 + *F3_n)*dt;
    *V4 = *V4 + 0.5*(*F4_0 + *F4_n)*dt;
    *V5 = *V5 + 0.5*(*F5_0 + *F5_n)*dt;
    *V6 = *V6 + 0.5*(*F6_0 + *F6_n)*dt;
    *V7 = *V7 + 0.5*(*F7_0 + *F7_n)*dt;
    *V8 = *V8 + 0.5*(*F8_0 + *F8_n)*dt;
    *V9 = *V9 + 0.5*(*F9_0 + *F9_n)*dt;

    *X1 = *X1_n, *X2 = *X2_n, *X3 = *X3_n, *X4 = *X4_n, *X5 = *X5_n, *X6 = *X6_n, *X7 = *X7_n, *X8 = *X8_n, *X9 = *X9_n;
    *F1_0 = *F1_n, *F2_0 = *F2_n, *F3_0 = *F3_n, *F4_0 = *F4_n, *F5_0 = *F5_n, *F6_0 = *F6_n, *F7_0 = *F7_n, *F8_0 = *F8_n, *F9_0 = *F9_n;
}

// Declare solution arrays
    matrix X1_sols[iterations/record],X2_sols[iterations/record],X3_sols[iterations/record],X4_sols[iterations/record];
    matrix X5_sols[iterations/record],X6_sols[iterations/record],X7_sols[iterations/record],X8_sols[iterations/record],X9_sols[iterations/record];
    matrix V1_sols[iterations/record],V2_sols[iterations/record],V3_sols[iterations/record],V4_sols[iterations/record];
    matrix V5_sols[iterations/record],V6_sols[iterations/record],V7_sols[iterations/record],V8_sols[iterations/record],V9_sols[iterations/record];
    matrix X1_sols_ds[iterations/record],X2_sols_ds[iterations/record],X3_sols_ds[iterations/record],X4_sols_ds[iterations/record];
    matrix X5_sols_ds[iterations/record],X6_sols_ds[iterations/record],X7_sols_ds[iterations/record],X8_sols_ds[iterations/record],X9_sols_ds[iterations/record];
    matrix V1_sols_ds[iterations/record],V2_sols_ds[iterations/record],V3_sols_ds[iterations/record],V4_sols_ds[iterations/record];
    matrix V5_sols_ds[iterations/record],V6_sols_ds[iterations/record],V7_sols_ds[iterations/record],V8_sols_ds[iterations/record],V9_sols_ds[iterations/record];

int main()
{
    // Declaring the matrices. g, and energy.
    long double g;
    long double Energy;
    matrix X1, X2, X3, X4, X5, X6, X7, X8, X9, X1_ds, X2_ds, X3_ds, X4_ds, X5_ds, X6_ds, X7_ds, X8_ds, X9_ds;
    matrix V1, V2, V3, V4, V5, V6, V7, V8, V9, V1_ds, V2_ds, V3_ds, V4_ds, V5_ds, V6_ds, V7_ds, V8_ds, V9_ds;

    // Loading Unperturbed and perturbed coordinates
    matrix X[18];
    std:: fstream thermalised_coordinates("Thermalised_branes_4x4[5740].txt", std:: ios:: in);
    thermalised_coordinates >> g;
    for (int i = 0; i < 18; i++)
    {
        for (int k = 0; k < N; k++)
        {
            for (int j = 0; j < N; j++)
            {
                thermalised_coordinates >> X[i](k, j);
            }
        }
    }
    thermalised_coordinates >> Energy;
    thermalised_coordinates.close();

    X1 = X[0], X2 = X[1], X3 = X[2], X4 = X[3], X5 = X[4], X6 = X[5], X7 = X[6], X8 = X[7], X9 = X[8];
    V1 = X[9], V2 = X[10], V3 = X[11], V4 = X[12], V5 = X[13], V6 = X[14], V7 = X[15], V8 = X[16], V9 = X[17];
    X1_ds = X1, X2_ds = X2, X3_ds = X3, X4_ds = X4, X5_ds = X5, X6_ds = X6, X7_ds = X7, X8_ds = X8, X9_ds = X9;
    V1_ds = V1, V2_ds = V2, V3_ds = V3, V4_ds = V4, V5_ds = V5, V6_ds = V6, V7_ds = V7, V8_ds = V8, V9_ds = V9;



    // Initializing F function at t = 0 for use in the update function
    matrix F1_0 = F(1,X1,X2,X3,X4,X5,X6,X7,X8,X9);
    matrix F2_0 = F(2,X1,X2,X3,X4,X5,X6,X7,X8,X9);
    matrix F3_0 = F(3,X1,X2,X3,X4,X5,X6,X7,X8,X9);
    matrix F4_0 = F(4,X1,X2,X3,X4,X5,X6,X7,X8,X9);
    matrix F5_0 = F(5,X1,X2,X3,X4,X5,X6,X7,X8,X9);
    matrix F6_0 = F(6,X1,X2,X3,X4,X5,X6,X7,X8,X9);
    matrix F7_0 = F(7,X1,X2,X3,X4,X5,X6,X7,X8,X9);
    matrix F8_0 = F(8,X1,X2,X3,X4,X5,X6,X7,X8,X9);
    matrix F9_0 = F(9,X1,X2,X3,X4,X5,X6,X7,X8,X9);

    matrix F1_0_ds = F(1,X1_ds, X2_ds, X3_ds, X4_ds, X5_ds, X6_ds, X7_ds, X8_ds, X9_ds);
    matrix F2_0_ds = F(2,X1_ds, X2_ds, X3_ds, X4_ds, X5_ds, X6_ds, X7_ds, X8_ds, X9_ds);
    matrix F3_0_ds = F(3,X1_ds, X2_ds, X3_ds, X4_ds, X5_ds, X6_ds, X7_ds, X8_ds, X9_ds);
    matrix F4_0_ds = F(4,X1_ds, X2_ds, X3_ds, X4_ds, X5_ds, X6_ds, X7_ds, X8_ds, X9_ds);
    matrix F5_0_ds = F(5,X1_ds, X2_ds, X3_ds, X4_ds, X5_ds, X6_ds, X7_ds, X8_ds, X9_ds);
    matrix F6_0_ds = F(6,X1_ds, X2_ds, X3_ds, X4_ds, X5_ds, X6_ds, X7_ds, X8_ds, X9_ds);
    matrix F7_0_ds = F(7,X1_ds, X2_ds, X3_ds, X4_ds, X5_ds, X6_ds, X7_ds, X8_ds, X9_ds);
    matrix F8_0_ds = F(8,X1_ds, X2_ds, X3_ds, X4_ds, X5_ds, X6_ds, X7_ds, X8_ds, X9_ds);
    matrix F9_0_ds = F(9,X1_ds, X2_ds, X3_ds, X4_ds, X5_ds, X6_ds, X7_ds, X8_ds, X9_ds);

    // Initializing X1_n's and F1_n's to pass addresses to the update function
    matrix X1_n, X2_n, X3_n, X4_n, X5_n, X6_n, X7_n, X8_n, X9_n, X1_n_ds, X2_n_ds, X3_n_ds, X4_n_ds, X5_n_ds, X6_n_ds, X7_n_ds, X8_n_ds, X9_n_ds;
    matrix F1_n, F2_n, F3_n, F4_n, F5_n, F6_n, F7_n, F8_n, F9_n, F1_n_ds, F2_n_ds, F3_n_ds, F4_n_ds, F5_n_ds, F6_n_ds, F7_n_ds, F8_n_ds, F9_n_ds;

    for (int i = 0; i < perturb_iterations; i++)
    {
        update(
            dt,
            &X1,&X2,&X3,&X4,&X5,&X6,&X7,&X8,&X9,
            &X1_n,&X2_n,&X3_n,&X4_n,&X5_n,&X6_n,&X7_n,&X8_n,&X9_n,
            &V1,&V2,&V3,&V4,&V5,&V6,&V7,&V8,&V9,
            &F1_0,&F2_0,&F3_0,&F4_0,&F5_0,&F6_0,&F7_0,&F8_0,&F9_0,
            &F1_n,&F2_n,&F3_n,&F4_n,&F5_n,&F6_n,&F7_n,&F8_n,&F9_n
        );

        perturb_update(
            dt,
            &X1_ds,&X2_ds,&X3_ds,&X4_ds,&X5_ds,&X6_ds,&X7_ds,&X8_ds,&X9_ds,
            &X1_n_ds,&X2_n_ds,&X3_n_ds,&X4_n_ds,&X5_n_ds,&X6_n_ds,&X7_n_ds,&X8_n_ds,&X9_n_ds,
            &V1_ds,&V2_ds,&V3_ds,&V4_ds,&V5_ds,&V6_ds,&V7_ds,&V8_ds,&V9_ds,
            &F1_0_ds,&F2_0_ds,&F3_0_ds,&F4_0_ds,&F5_0_ds,&F6_0_ds,&F7_0_ds,&F8_0_ds,&F9_0_ds,
            &F1_n_ds,&F2_n_ds,&F3_n_ds,&F4_n_ds,&F5_n_ds,&F6_n_ds,&F7_n_ds,&F8_n_ds,&F9_n_ds
        );

        if (i % record == 0)
        {
            X1_sols[i/record] = X1, X2_sols[i/record] = X2, X3_sols[i/record] = X3, X4_sols[i/record] = X4, 
            X5_sols[i/record] = X5, X6_sols[i/record] = X6, X7_sols[i/record] = X7, X8_sols[i/record] = X8, X9_sols[i/record] = X9;
            
            X1_sols_ds[i/record] = X1_ds, X2_sols_ds[i/record] = X2_ds, X3_sols_ds[i/record] = X3_ds, X4_sols_ds[i/record] = X4_ds,
            X5_sols_ds[i/record] = X5_ds, X6_sols_ds[i/record] = X6_ds, X7_sols_ds[i/record] = X7_ds, X8_sols_ds[i/record] = X8_ds,
            X9_sols_ds[i/record] = X9_ds;

            V1_sols[i/record] = V1, V2_sols[i/record] = V2, V3_sols[i/record] = V3, V4_sols[i/record] = V4,
            V5_sols[i/record] = V5, V6_sols[i/record] = V6, V7_sols[i/record] = V7, V8_sols[i/record] = V8,
            V9_sols[i/record] = V9;

            V1_sols_ds[i/record] = V1_ds, V2_sols_ds[i/record] = V2_ds, V3_sols_ds[i/record] = V3_ds, V4_sols_ds[i/record] = V4_ds,
            V5_sols_ds[i/record] = V5_ds, V6_sols_ds[i/record] = V6_ds, V7_sols_ds[i/record] = V7_ds, V8_sols_ds[i/record] = V8_ds,
            V9_sols_ds[i/record] = V9_ds;

        }
/*
        if (i%1000 == 0)
        {
            std:: cout << std::setprecision(15) << i << std::endl;
            std:: cout << "time: " << std::time(nullptr) - start << std:: endl;
            std:: cout << H(g,X1,X2,X3,X4,X5,X6,X7,X8,X9,V1,V2,V3,V4,V5,V6,V7,V8,V9) << std:: endl;
            std:: cout << H(g,X1_ds,X2_ds,X3_ds,X4_ds,X5_ds,X6_ds,X7_ds,X8_ds,X9_ds,V1_ds,V2_ds,V3_ds,V4_ds,V5_ds,V6_ds,V7_ds,V8_ds,V9_ds) << std:: endl;
        }
*/
    }
            std:: cout <<std::setprecision(15)<< "Ideal perturb finished" << '\n';
            std:: cout << "X1: " << X1 << std:: endl << "X1_ds: " << X1_ds << std:: endl << "X2: " << X2 << std:: endl << "X2_ds: " << X2_ds <<std:: endl;
            std:: cout << "X3: " << X3 << std:: endl << "X3_ds: " << X3_ds << std:: endl << "X4: " << X4 << std:: endl << "X4_ds: " << X4_ds <<std:: endl;
            std:: cout << "X5: " << X5 << std:: endl << "X5_ds: " << X5_ds << std:: endl << "X6: " << X6 << std:: endl << "X6_ds: " << X6_ds <<std:: endl;
            std:: cout << "X7: " << X7 << std:: endl << "X7_ds: " << X7_ds << std:: endl << "X8: " << X8 << std:: endl << "X8_ds: " << X8_ds <<std:: endl;
            std:: cout << "X9: " << X9 << std:: endl << "X9_ds: " << X9_ds << std:: endl;
      

    // Run update function
    for (int i = 0; i < 600; i++)
    {
        update(
            dt,
            &X1,&X2,&X3,&X4,&X5,&X6,&X7,&X8,&X9,
            &X1_n,&X2_n,&X3_n,&X4_n,&X5_n,&X6_n,&X7_n,&X8_n,&X9_n,
            &V1,&V2,&V3,&V4,&V5,&V6,&V7,&V8,&V9,
            &F1_0,&F2_0,&F3_0,&F4_0,&F5_0,&F6_0,&F7_0,&F8_0,&F9_0,
            &F1_n,&F2_n,&F3_n,&F4_n,&F5_n,&F6_n,&F7_n,&F8_n,&F9_n);
    
        update(
            dt,
            &X1_ds,&X2_ds,&X3_ds,&X4_ds,&X5_ds,&X6_ds,&X7_ds,&X8_ds,&X9_ds,
            &X1_n_ds,&X2_n_ds,&X3_n_ds,&X4_n_ds,&X5_n_ds,&X6_n_ds,&X7_n_ds,&X8_n_ds,&X9_n_ds,
            &V1_ds,&V2_ds,&V3_ds,&V4_ds,&V5_ds,&V6_ds,&V7_ds,&V8_ds,&V9_ds,
            &F1_0_ds,&F2_0_ds,&F3_0_ds,&F4_0_ds,&F5_0_ds,&F6_0_ds,&F7_0_ds,&F8_0_ds,&F9_0_ds,
            &F1_n_ds,&F2_n_ds,&F3_n_ds,&F4_n_ds,&F5_n_ds,&F6_n_ds,&F7_n_ds,&F8_n_ds,&F9_n_ds);

        if (i % record == 0)
        {
            
            X1_sols[i/record] = X1, X2_sols[i/record] = X2, X3_sols[i/record] = X3, X4_sols[i/record] = X4, 
            X5_sols[i/record] = X5, X6_sols[i/record] = X6, X7_sols[i/record] = X7, X8_sols[i/record] = X8, X9_sols[i/record] = X9;
            
            X1_sols_ds[i/record] = X1_ds, X2_sols_ds[i/record] = X2_ds, X3_sols_ds[i/record] = X3_ds, X4_sols_ds[i/record] = X4_ds,
            X5_sols_ds[i/record] = X5_ds, X6_sols_ds[i/record] = X6_ds, X7_sols_ds[i/record] = X7_ds, X8_sols_ds[i/record] = X8_ds,
            X9_sols_ds[i/record] = X9_ds;

            V1_sols[i/record] = V1, V2_sols[i/record] = V2, V3_sols[i/record] = V3, V4_sols[i/record] = V4,
            V5_sols[i/record] = V5, V6_sols[i/record] = V6, V7_sols[i/record] = V7, V8_sols[i/record] = V8,
            V9_sols[i/record] = V9;

            V1_sols_ds[i/record] = V1_ds, V2_sols_ds[i/record] = V2_ds, V3_sols_ds[i/record] = V3_ds, V4_sols_ds[i/record] = V4_ds,
            V5_sols_ds[i/record] = V5_ds, V6_sols_ds[i/record] = V6_ds, V7_sols_ds[i/record] = V7_ds, V8_sols_ds[i/record] = V8_ds,
            V9_sols_ds[i/record] = V9_ds;
        }

        if (i%200 == 0)
        {
            std:: cout << i << std::endl;
            std:: cout << H(g,X1,X2,X3,X4,X5,X6,X7,X8,X9,V1,V2,V3,V4,V5,V6,V7,V8,V9) << std:: endl;
            std:: cout << H(g,X1_ds,X2_ds,X3_ds,X4_ds,X5_ds,X6_ds,X7_ds,X8_ds,X9_ds,V1_ds,V2_ds,V3_ds,V4_ds,V5_ds,V6_ds,V7_ds,V8_ds,V9_ds) << std:: endl;
    
            std:: cout <<std::setprecision(15)<< "X1: " << X1 << std:: endl << "X1_ds: " << X1_ds << std:: endl << "X2: " << X2 << std:: endl << "X2_ds: " << X2_ds <<std:: endl;
            std:: cout << "X3: " << X3 << std:: endl << "X3_ds: " << X3_ds << std:: endl << "X4: " << X4 << std:: endl << "X4_ds: " << X4_ds <<std:: endl;
            std:: cout << "X5: " << X5 << std:: endl << "X5_ds: " << X5_ds << std:: endl << "X6: " << X6 << std:: endl << "X6_ds: " << X6_ds <<std:: endl;
            std:: cout << "X7: " << X7 << std:: endl << "X7_ds: " << X7_ds << std:: endl << "X8: " << X8 << std:: endl << "X8_ds: " << X8_ds <<std:: endl;
            std:: cout << "X9: " << X9 << std:: endl << "X9_ds: " << X9_ds << std:: endl;
            
    
        }
    }
    
    // Writing to files
    std:: fstream X_file("C:/Users/cilli/DIAS.cpp/D0_Brane_Solutions/6x6/E_183_g/X_sols.txt",std:: ios:: out);
    std:: fstream X_file_ds("C:/Users/cilli/DIAS.cpp/D0_Brane_Solutions/6x6/E_183_g/X_sols_ds.txt", std:: ios:: out);
    std:: fstream V_file("C:/Users/cilli/DIAS.cpp/D0_Brane_Solutions/6x6/E_183_g/V_sols.txt", std:: ios:: out);
    std:: fstream V_file_ds("C:/Users/cilli/DIAS.cpp/D0_Brane_Solutions/6x6/E_183_g/V_sols_ds.txt", std:: ios:: out);

    for (int k = 0; k < iterations/record; k ++)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    X_file << X1_sols[k](i,j);
                    X_file_ds << X1_sols_ds[k](i,j);
                    V_file << V1_sols[k](i,j);
                    V_file_ds << V1_sols_ds[k](i,j);
                }
            }
            X_file << ":";
            X_file_ds << ":";
            V_file << ":";
            V_file_ds << ":";
        }
    for (int k = 0; k < iterations/record; k ++)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    X_file << X2_sols[k](i,j);
                    X_file_ds << X2_sols_ds[k](i,j);
                    V_file << V2_sols[k](i,j);
                    V_file_ds << V2_sols_ds[k](i,j);
                }
            }
            X_file << ":";
            X_file_ds << ":";
            V_file << ":";
            V_file_ds << ":";
        }
    for (int k = 0; k < iterations/record; k ++)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    X_file << X3_sols[k](i,j);
                    X_file_ds << X3_sols_ds[k](i,j);
                    V_file << V3_sols[k](i,j);
                    V_file_ds << V3_sols_ds[k](i,j);
                }
            }
            X_file << ":";
            X_file_ds << ":";
            V_file << ":";
            V_file_ds << ":";
        }
    for (int k = 0; k < iterations/record; k ++)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    X_file << X4_sols[k](i,j);
                    X_file_ds << X4_sols_ds[k](i,j);
                    V_file << V4_sols[k](i,j);
                    V_file_ds << V4_sols_ds[k](i,j);
                }
            }
            X_file << ":";
            X_file_ds << ":";
            V_file << ":";
            V_file_ds << ":";
        }
    for (int k = 0; k < iterations/record; k ++)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    X_file << X5_sols[k](i,j);
                    X_file_ds << X5_sols_ds[k](i,j);
                    V_file << V5_sols[k](i,j);
                    V_file_ds << V5_sols_ds[k](i,j);
                }
            }
            X_file << ":";
            X_file_ds << ":";
            V_file << ":";
            V_file_ds << ":";
        }
    for (int k = 0; k < iterations/record; k ++)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    X_file << X6_sols[k](i,j);
                    X_file_ds << X6_sols_ds[k](i,j);
                    V_file << V6_sols[k](i,j);
                    V_file_ds << V6_sols_ds[k](i,j);
                }
            }
            X_file << ":";
            X_file_ds << ":";
            V_file << ":";
            V_file_ds << ":";
        }
    for (int k = 0; k < iterations/record; k ++)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    X_file << X7_sols[k](i,j);
                    X_file_ds << X7_sols_ds[k](i,j);
                    V_file << V7_sols[k](i,j);
                    V_file_ds << V7_sols_ds[k](i,j);
                }
            }
            X_file << ":";
            X_file_ds << ":";
            V_file << ":";
            V_file_ds << ":";
        }
    for (int k = 0; k < iterations/record; k ++)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    X_file << X8_sols[k](i,j);
                    X_file_ds << X8_sols_ds[k](i,j);
                    V_file << V8_sols[k](i,j);
                    V_file_ds << V8_sols_ds[k](i,j);
                }
            }
            X_file << ":";
            X_file_ds << ":";
            V_file << ":";
            V_file_ds << ":";
        }
    for (int k = 0; k < iterations/record; k ++)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    X_file << X9_sols[k](i,j);
                    X_file_ds << X9_sols_ds[k](i,j);
                    V_file << V9_sols[k](i,j);
                    V_file_ds << V9_sols_ds[k](i,j);
                }
            }
            X_file << ":";
            X_file_ds << ":";
            V_file << ":";
            V_file_ds << ":";
        }
    
    X_file.close();
    X_file_ds.close();
    V_file.close();
    V_file_ds.close();
    return 0;
}