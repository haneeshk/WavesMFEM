#pragma once
#include <string>
#include "mfem.hpp"
using namespace mfem; 

void saveData(std::string, const int, const Vector&, const GridFunction &);