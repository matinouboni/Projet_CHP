#ifndef _FUNCTION_CPP

#include "Function.h"
#include <cmath>

Function::Function(DataFile* data_file) :
_df(data_file)
{
   
}

double Function::Initial_condition(const double x, const double y) const
{
   if (this->_df->Get_cas() == 0)
   {
      double pi(3.14159265358979323846264338327950288);
      return sin(2.0*pi*x/(_df->Get_xmax()-_df->Get_xmin()))*sin(2.0*pi*y/(_df->Get_ymax()-_df->Get_ymin()));
   }
   else if (this->_df->Get_cas() == 1)
   {
      return exp(-(x*x)/(0.0075)-(y*y)/(0.0075));
   }
   else if (this->_df->Get_cas() == 2)
   {  
      if (sqrt((x-0.0)*(x-0.0)+(y-0.0)*(y-0.0)) <= 0.4){ 
         return 1.0;
      }
      else{
         return 0.0;
      }
   }
   else if (this->_df->Get_cas() == 3)
   {  
      if (sqrt((x-0.0)*(x-0.0)+(y-0.0)*(y-0.0)) <= 0.4){ 
         return 1.0;
      }
      else{
         return 0.0;
      }
   }
   else if (this->_df->Get_cas() == 4)
   {  
      if ((sqrt((x-0.5)*(x-0.5)+(y-0.75)*(y-0.75)) <= 0.15) && ((0.025 <= fabs(x-0.5)) || (0.85 <= y))){ 
         return 1.0;
      }
      else{
         return 0.0;
      }
   }
   else if (this->_df->Get_cas() == 5)
   {  
      if (sqrt((x-0.5)*(x-0.5)+(y-0.75)*(y-0.75)) <= 0.15){ 
         return 1.0;
      }
      else{
         return 0.0;
      }
   }
   else
   {
      return 0.0;
   }
}

double Function::Velocity_x(const double x, const double y, const double t) const
{
   if (this->_df->Get_cas() == 0)
   {
      return 1.0;
   }
   else if (this->_df->Get_cas() == 1)
   {
      return 1.0;
   }
   else if (this->_df->Get_cas() == 2)
   {
      return 0.5;
   }
   else if (this->_df->Get_cas() == 3)
   {
      return -y;
   }
   else if (this->_df->Get_cas() == 4)
   {
      return 0.5-y;
   }
   else if (this->_df->Get_cas() == 5)
   {
      double pi(3.14159265358979323846264338327950288);
      return 2.0*sin(pi*x)*sin(pi*x)*sin(2.0*pi*y)*cos(2.0*pi*(t/6.0));
   }
   else
   {
      return 0.0;
   }
   
}


double Function::Velocity_y(const double x, const double y, const double t) const
{
   if (this->_df->Get_cas() == 0)
   {
      return 0.0;
   }
   else if (this->_df->Get_cas() == 1)
   {
      return 0.0;
   }
   else if (this->_df->Get_cas() == 2)
   {
      return 0.5;
   }
   else if (this->_df->Get_cas() == 3)
   {
      return x;
   }
   else if (this->_df->Get_cas() == 4)
   {
      return x-0.5;
   }
   else if (this->_df->Get_cas() == 5)
   {
      double pi(3.14159265358979323846264338327950288);
      return -sin(2.0*pi*x)*sin(pi*y)*sin(pi*y)*cos(2.0*pi*(t/6.0));
   }
   else
   {
      return 0.0;
   }
}

double Function::Exact_solution(const double x, const double y, const double t) const
{
  return Initial_condition(x - Velocity_x(x,y,t)*t, y - Velocity_y(x,y,t)*t);
}


#define _FUNCTION_CPP
#endif
