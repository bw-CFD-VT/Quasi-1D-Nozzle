// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#ifndef TIMESTEP_HPP
#define TIMESTEP_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window

using namespace std;

//---------------------------------------- Time Step, dt ----------------------------------------------------------------//
void Time_Step (int imax, double CFL, double dx, vector<vector<vector<double> > > V_cell_center,
                vector<double> &lambda_max, vector<double> &a,vector<double> &dt);
//-----------------------------------------------------------------------------------------------------------------------//

#endif