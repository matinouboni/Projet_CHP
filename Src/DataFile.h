#ifndef _DATA_FILE_H

#include <string>
#include <vector>
#include <iostream>
// Définition de la classe

class DataFile {

   private:
      const std::string _file_name;
      int _cas, _Nx, _Ny, _space_scheme, _time_scheme; 
      double _xmin, _xmax, _ymin, _ymax, _Tf, _CFL;
      double _dx, _dy, _dt;


   public: // Méthodes et opérateurs de la classe
   DataFile(std::string file_name);

   const int Get_cas() const {return _cas;};
   const int Get_Nx() const {return _Nx;};
   const int Get_Ny() const {return _Ny;};
   const int Get_space_scheme() const {return _space_scheme;};
   const int Get_time_scheme() const {return _time_scheme;};
   const double Get_xmin() const {return _xmin;};
   const double Get_xmax() const {return _xmax;};
   const double Get_ymin() const {return _ymin;};
   const double Get_ymax() const {return _ymax;};
   const double Get_Tf() const {return _Tf;};
   const double Get_CFL() const {return _CFL;};
   const double Get_dx() const {return (_xmax - _xmin)/(double(_Nx));};
   const double Get_dy() const {return (_ymax - _ymin)/(double(_Ny));};
   const double Get_dt() const {return _dt;};
   
};

#define _DATA_FILE_H
#endif
