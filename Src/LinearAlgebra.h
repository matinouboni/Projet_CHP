#ifndef _LINEAR_ALGEBRA_H

#include <vector>
#include "mpi.h"

double DotProduct(const std::vector<double> a, const std::vector<double> b, const int &Nx, const int &Ny, const int &Me, const int &Np, MPI_Status &status);
std::vector<double> MultiplyBy(const std::vector<double> a, const double b);
std::vector<double> AddVector(const std::vector<double> a, const std::vector<double> b);
std::vector<double> SubVector(const std::vector<double> a, const std::vector<double> b);
void charge(int me, int n, int np, int *iBeg, int *iEnd);

#define _LINEAR_ALGEBRA_H
#endif
