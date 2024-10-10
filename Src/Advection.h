#ifndef _ADVECTION_H

#include <string>
#include "Function.h"
#include <vector>
#include <cmath>
#include "mpi.h"

class Advection
{
private:
	Function* _fct;
	DataFile* _df;
	double _dt;
public:
	// Constructeur
	Advection(Function* function, DataFile* data_file);
	std::vector<double> MatVecProd(const std::vector<double> &U, const double t, const int &iBeg, const int &iEnd, const int &Me, const int &Np, MPI_Status &status);
	void InitialCondition(std::vector<double> &U, const int &iBeg, const int &iEnd);
	const double Get_dt() const {return _dt;};
	std::vector<double> ExactSol(const double t, const int &iBeg, const int &iEnd);
};

#define _ADVECTION_H
#endif
