#ifndef _ADVECTION_CPP

#include "Advection.h"
#include <iostream>
#include<cmath>

using namespace std;

// Constructeur
Advection::Advection(Function* function, DataFile* data_file) :
_fct(function), _df(data_file)
{
    // Definition du pas de temps si le schema est explicite decentre
    if (_df->Get_time_scheme() == 1){
        // Cas test 0, 1 ou 2
        if((_df->Get_cas() == 1) || (_df->Get_cas() == 2) || (_df->Get_cas() == 0)){
            _dt = _df->Get_CFL()/(fabs(_fct->Velocity_x(0,0,0))/_df->Get_dx()+fabs(_fct->Velocity_y(0,0,0))/_df->Get_dy());
        }
        // Cas test 3
        else if (_df->Get_cas() == 3){
            _dt = _df->Get_CFL()/(1.0/_df->Get_dx()+1.0/_df->Get_dy());
        }
        // Cas test 4
        else if (_df->Get_cas() == 4){
            _dt = _df->Get_CFL()/(0.5/_df->Get_dx()+0.5/_df->Get_dy());
        }
        // Cas test 5
        else if (_df->Get_cas() == 5){
            _dt = _df->Get_CFL()/(2.0/_df->Get_dx()+1.0/_df->Get_dy());
        }
    }
    // Definition du pas de temps si le schema est implite centre
    else if (_df->Get_time_scheme() == 2){
        // Si cas test 0 pour retrouver l'ordre 2 du schema centre
        if (_df->Get_cas() != 0){
             _dt = _df->Get_dt();
        }
        else{
             _dt = (1.0/(fabs(_fct->Velocity_x(0,0,0))/_df->Get_dx()+fabs(_fct->Velocity_y(0,0,0))/_df->Get_dy()))*(1.0/(fabs(_fct->Velocity_x(0,0,0))/_df->Get_dx()+fabs(_fct->Velocity_y(0,0,0))/_df->Get_dy()));
        }
    }
}

void Advection::InitialCondition(std::vector<double> &U, const int &iBeg, const int &iEnd){
    // Calcul du vecteur solution initiale
    int i, j;
    for (int k = 0; k < iEnd-iBeg+1; ++k){
        i = (k+iBeg)%_df->Get_Nx() + 1;
        j = (k+iBeg)/_df->Get_Nx() + 1;
        U[k] = _fct->Initial_condition(_df->Get_xmin() + i*_df->Get_dx(),_df->Get_ymin() + j*_df->Get_dy());
    }
}

std::vector<double> Advection::MatVecProd(const std::vector<double> &U, const double t, const int &iBeg, const int &iEnd, const int &Me, const int &Np, MPI_Status &status){
    // Produit matrice vecteur creux
    int i, j;
    std::vector<double> X(U.size());
    int Nx(_df->Get_Nx()), Ny(_df->Get_Ny());
    double dx(_df->Get_dx()), dy(_df->Get_dy()), xmin(_df->Get_xmin()), ymin(_df->Get_ymin()); 
    int Tag11(11), Tag12(12);
    std::vector<double> U_Mem1(Nx,0.0), U_Mep1(Nx,0.0);
    // Envoie et reception des portions de veceteur des proc Me-1 et Me+1 par le proc Me
    MPI_Send(&U[iEnd-iBeg-Nx+1], Nx, MPI_DOUBLE, (Me+1)%Np, Tag12, MPI_COMM_WORLD);
    MPI_Recv(&U_Mem1[0], Nx, MPI_DOUBLE, (Me+Np-1)%Np, Tag12, MPI_COMM_WORLD, &status);
    MPI_Send(&U[0], Nx, MPI_DOUBLE, (Me+Np-1)%Np, Tag11, MPI_COMM_WORLD);
    MPI_Recv(&U_Mep1[0], Nx, MPI_DOUBLE, (Me+1)%Np, Tag11, MPI_COMM_WORLD, &status);
    std::vector<double> U_extend(U.size(),0.0);
    // Reconstrcution d'un vecteur prennant en compte les portions de vecteur de Me-1 et Me+1
    U_extend = U;
    U_extend.insert(U_extend.end(), U_Mep1.begin(), U_Mep1.end());
    U_extend.insert(U_extend.begin(), U_Mem1.begin(), U_Mem1.end());

    for (int k = Nx; k < U_extend.size()-Nx; ++k){
        i = (k+iBeg-Nx)%Nx + 1;
        j = (k+iBeg-Nx)/Nx + 1;
        double Vx(_fct->Velocity_x(xmin + i*dx, ymin + j*dy, t));
        double Vy(_fct->Velocity_y(xmin + i*dx, ymin + j*dy, t));
        double Vx_p((Vx+fabs(Vx))/2.0), Vx_m((Vx-fabs(Vx))/2.0), Vy_p((Vy+fabs(Vy))/2.0), Vy_m((Vy-fabs(Vy))/2.0);
        // Si schema explcite decentre
        if ((_df->Get_time_scheme() == 1) && (_df->Get_space_scheme() == 1)){
            double alpha(1.0 - _dt/dx*(Vx_p-Vx_m) - _dt/dy*(Vy_p-Vy_m));
            double beta(_dt*Vx_p/dx);
            double b(-_dt*Vx_m/dx);
            double c(-_dt*Vy_m/dy);
            double gamma(_dt*Vy_p/dy);

            if (i == 1){
                X[k-Nx] = gamma*U_extend[k-Nx] + beta*U_extend[k-1+Nx] + alpha*U_extend[k] + b*U_extend[k+1] + c*U_extend[k+Nx];
            }
            else if (i == Nx){    
                X[k-Nx] = gamma*U_extend[k-Nx] + beta*U_extend[k-1] + alpha*U_extend[k] + b*U_extend[k+1-Nx] + c*U_extend[k+Nx];
            }
            else {
                X[k-Nx] = gamma*U_extend[k-Nx] + beta*U_extend[k-1] + alpha*U_extend[k] + b*U_extend[k+1] + c*U_extend[k+Nx];
            }

        }
        // Si schema implicite centre
        else if ((_df->Get_time_scheme() == 2) && (_df->Get_space_scheme() == 2)){
                double alpha(1.0);
                double beta(-_dt*Vx/(2.0*dx));
                double b(_dt*Vx/(2.0*dx));
                double c(_dt*Vy/(2.0*dy));
                double gamma(-_dt*Vy/(2.0*dy));

                if (i == 1){
                    X[k-Nx] = gamma*U_extend[k-Nx] + beta*U_extend[k-1+Nx] + alpha*U_extend[k] + b*U_extend[k+1] + c*U_extend[k+Nx];
                }
                else if (i == Nx){    
                    X[k-Nx] = gamma*U_extend[k-Nx] + beta*U_extend[k-1] + alpha*U_extend[k] + b*U_extend[k+1-Nx] + c*U_extend[k+Nx];
                }
                else {
                    X[k-Nx] = gamma*U_extend[k-Nx] + beta*U_extend[k-1] + alpha*U_extend[k] + b*U_extend[k+1] + c*U_extend[k+Nx];
                }
        }
    }
    return X;
}

std::vector<double> Advection::ExactSol(const double t, const int &iBeg, const int &iEnd){
    // Calcul de la solution exacte dans un vecteur
    std::vector<double> ExactSol(iEnd-iBeg+1);
    int i, j;
        for (int k = 0; k < iEnd-iBeg+1; ++k){
            i = (k+iBeg)%_df->Get_Nx() + 1;
            j = (k+iBeg)/_df->Get_Nx() + 1;
            ExactSol[k] = _fct->Exact_solution(_df->Get_xmin() + i*_df->Get_dx(),_df->Get_ymin() + j*_df->Get_dy(),t);
        }
        return ExactSol;
}

#define _ADVECTION_CPP
#endif
