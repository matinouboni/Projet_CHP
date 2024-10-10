#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"

TimeScheme::TimeScheme(DataFile* data_file, Advection* adv) :
_df(data_file), _adv(adv)
{

}

void TimeScheme::SaveSol(const std::vector<double> &U, std::string n_sol, int n, const int &iBeg, const int &iEnd, const  int &Me){
    // Sauvegarde de la solution
    std::string n_file = "../Results/" + n_sol + "." + std::to_string(n) + "." + "00" + std::to_string(Me) +".dat";
    // std::ofstream monflux;
    // monflux.open(n_file, std::ios::out);
    FILE* monflux;
    const char *cstr = n_file.c_str();
    monflux = fopen(cstr, "w");
    int i, j;
    for (int k = 0; k < iEnd-iBeg+1; ++k){
        i = (k+iBeg)%_df->Get_Nx() + 1;
        j = (k+iBeg)/_df->Get_Nx()+ 1;
            // << est plus lent que fprintf
            // monflux << _df->Get_xmin() + i*_df->Get_dx() << " " << _df->Get_ymin() + j*_df->Get_dy() << " " << U[k] << std::endl; 
            fprintf(monflux, "%lf   %lf   %lf\n", _df->Get_xmin() + i*_df->Get_dx(), _df->Get_ymin() + j*_df->Get_dy(), U[k]);
    }
    // monflux.close();
    fclose(monflux);
}

ExplicitScheme::ExplicitScheme(DataFile* data_file, Advection* adv) : 
TimeScheme(data_file, adv)
{

}

void ExplicitScheme::Integrate(double &t, std::vector<double> &U, const int &iBeg, const int &iEnd, const int &Me, const int &Np, MPI_Status &status){
    // Calcul de la solution a l'instant n+1
    U = _adv->MatVecProd(U, t, iBeg, iEnd, Me, Np, status);
    t += _adv->Get_dt();
}


ImplicitScheme::ImplicitScheme(DataFile* data_file, Advection* adv) : 
TimeScheme(data_file, adv)
{

}

std::vector<double> ImplicitScheme::BiCGstab(std::vector<double> &U, double &t, const int &iBeg, const int &iEnd, const int &Me, const int &Np, MPI_Status &status){
    // Algorithme du BiCG stab
    std::vector<double> r(U.size()), r_old(U.size()), r_tilde(U.size()), p(U.size()), nu(U.size()), h(U.size()), s(U.size());
    std::vector<double> t1(U.size()), x(U.size()), p_old(U.size());
    double rho(0.0), rho_old(0.0), alpha(0.0), omega(0.0), beta(0.0);
    int Nmax(10000), k(0);
    int Nx(_df->Get_Nx()), Ny(_df->Get_Ny());
    x = U;
    r = SubVector(U,_adv->MatVecProd(x, t, iBeg, iEnd, Me, Np, status));
    r_tilde = r;
    rho = DotProduct(r_tilde,r,Nx,Ny,Me,Np,status);
    p = r;
    while (k <= Nmax){
        nu = _adv->MatVecProd(p, t, iBeg, iEnd, Me, Np, status);
        alpha = rho/DotProduct(r_tilde,nu,Nx,Ny,Me,Np,status);
        h = AddVector(x,MultiplyBy(p,alpha));
        s = SubVector(r,MultiplyBy(nu,alpha));
        if (sqrt(DotProduct(s,s,Nx,Ny,Me,Np,status)) <= 1.e-12){
            x = h;
            return x;
        }
        t1 = _adv->MatVecProd(s, t, iBeg, iEnd, Me, Np, status);
        omega = DotProduct(t1,s,Nx,Ny,Me,Np,status)/DotProduct(t1,t1,Nx,Ny,Me,Np,status);
        x = AddVector(h,MultiplyBy(s,omega));
        r_old = r;
        r = SubVector(s,MultiplyBy(t1,omega));
        if (sqrt(DotProduct(r,r,Nx,Ny,Me,Np,status)) <= 1.e-12) {
            return x;
        }
        rho_old = rho;
        rho = DotProduct(r_tilde,r,Nx,Ny,Me,Np,status);
        beta = (rho/rho_old)*(alpha/omega);
        p_old = p;
        p = AddVector(r,MultiplyBy(SubVector(p_old,MultiplyBy(nu,omega)),beta));
        k++;
    }
    if (k >= Nmax){
        printf("Pas de convergence");
        exit(0);
    }
    return x;
}

void ImplicitScheme::Integrate(double &t, std::vector<double> &U, const int &iBeg, const int &iEnd, const int &Me, const int &Np, MPI_Status &status){
    // Calcul de la solution a l'instant n+1
    t += _adv->Get_dt();
    U = BiCGstab(U, t, iBeg, iEnd, Me, Np, status);
}

#define _TIME_SCHEME_CPP
#endif