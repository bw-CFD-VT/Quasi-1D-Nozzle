// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#ifndef FLUX_HPP
#define FLUX_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window

using namespace std;


//---------------------------------------- Flux Vector, F @ iteration = counter -----------------------------------------//
void Flux (int counter,int NI, vector<vector<vector<double> > > V_Boundary, vector<vector<vector<double> > > U_cell_center,
           vector<vector<vector<double> > > &F);
//-----------------------------------------------------------------------------------------------------------------------//


#endif
