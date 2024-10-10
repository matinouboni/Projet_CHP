#ifndef _DATA_FILE_CPP

#include "DataFile.h"
#include "../Data/toml.hpp"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

DataFile::DataFile(std::string file_name)
: _file_name(file_name)
{
   // Lecture du fichier de données
   auto config = toml::parse(file_name);

   // Other
   const auto& parameter = toml::find(config, "parameter");
   this->_cas = toml::find<int>(parameter, "cas");
   this->_xmin = toml::find<double>(parameter, "xmin");
   this->_xmax = toml::find<double>(parameter, "xmax");
   this->_ymin = toml::find<double>(parameter, "ymin");
   this->_ymax = toml::find<double>(parameter, "ymax");
   this->_Tf = toml::find<double>(parameter, "Tf");
   this->_Nx = toml::find<int>(parameter, "Nx");
   this->_Ny = toml::find<int>(parameter, "Ny");
   this->_CFL = toml::find<double>(parameter, "CFL");
   this->_space_scheme = toml::find<int>(parameter, "space_scheme");
   this->_time_scheme = toml::find<int>(parameter, "time_scheme");
   this->_dt = toml::find<double>(parameter, "dt");

}

#define _DATA_FILE_CPP
#endif
