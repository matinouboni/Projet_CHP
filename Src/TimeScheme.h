#ifndef _TIME_SCHEME_H

#include "Advection.h"
#include "LinearAlgebra.h"
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>


class TimeScheme
{
    protected:
    DataFile* _df;
    Advection* _adv;

    public:
    TimeScheme(DataFile* data_file, Advection* adv);
    virtual void Integrate(double &t, std::vector<double> &U, const int &iBeg, const int &iEnd, const int &Me, const int &Np, MPI_Status &status) = 0;
    void SaveSol(const std::vector<double> &U, std::string n_sol, int n, const int &iBeg, const int &iEnd, const  int &Me);
};

class ExplicitScheme : public TimeScheme
{
    public:
    ExplicitScheme(DataFile* data_file, Advection* adv);
    void Integrate(double &t, std::vector<double> &U, const int &iBeg, const int &iEnd, const int &Me, const int &Np, MPI_Status &status);
};

class ImplicitScheme : public TimeScheme
{
    public:
    ImplicitScheme(DataFile* data_file, Advection* adv);
    void Integrate(double &t, std::vector<double> &U, const int &iBeg, const int &iEnd, const int &Me, const int &Np, MPI_Status &status);
    std::vector<double> BiCGstab(std::vector<double> &U, double &t, const int &iBeg, const int &iEnd, const int &Me, const int &Np, MPI_Status &status);
};

#define _TIME_SCHEME_H
#endif