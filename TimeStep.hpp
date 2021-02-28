#ifndef TIMESTEP_HPP
#define TIMESTEP_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window

using namespace std;

void time_step (double CFL, double dx, vector<vector<vector<double>>> V_cell_center, vector<double> &dt);

#endif