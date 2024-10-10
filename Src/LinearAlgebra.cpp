#ifndef _LINEAR_ALGEBRA_CPP

#include "LinearAlgebra.h"

double DotProduct(const std::vector<double> a, const std::vector<double> b, const int &Nx, const int &Ny, const int &Me, const int &Np, MPI_Status &status) {
    // Produit scalaire
    double result(0.0), result_loc(0.0);
    int Tag1(11);
    int iBeg, iEnd, jBeg, jEnd;
    std::vector<double> A(Nx*Ny), B(Nx*Ny);
    
    charge(Me, Nx*Ny, Np, &iBeg, &iEnd);

    for (int i = 0; i < iEnd-iBeg+1; i++){
        result_loc += a[i]*b[i];
    }

    // Reduction
    MPI_Allreduce(&result_loc, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return result;
}


std::vector<double> MultiplyBy(const std::vector<double> a, const double b){
    // Multiplication d'un vecteur a par un scalaire b
    std::vector<double> result(a.size());
    for(int i = 0; i < a.size(); i++){
        result[i] = b*a[i];
    }
    return result;
}

std::vector<double> AddVector(const std::vector<double> a, const std::vector<double> b){
    // Somme de deux vecteurs
    std::vector<double> result(a.size());
    for(int i = 0; i < a.size(); i++){
        result[i] = a[i] + b[i];
    }
    return result;
}

std::vector<double> SubVector(const std::vector<double> a, const std::vector<double> b){
    // Difference de deux vecteurs
    std::vector<double> result(a.size());
    for(int i = 0; i < a.size(); i++){
        result[i] = a[i] - b[i];
    }
    return result;
}

void charge(int me, int n, int np, int *iBeg, int *iEnd){
    // Fonction charge
    if (me < n%np){
        *iBeg = (n/np+1)*me;
        *iEnd = (n/np+1)*(me+1) - 1;
    }else{
        *iBeg = n%np*(n/np+1) + (me - n%np)*(n/np);
        *iEnd = *iBeg + n/np - 1;
    }
}

#define _LINEAR_ALGEBRA_CPP
#endif