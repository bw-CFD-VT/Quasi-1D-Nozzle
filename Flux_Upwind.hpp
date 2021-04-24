// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#ifndef FLUX_UPWIND_HPP
#define FLUX_UPWIND_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window

using namespace std;

//---------------------------------------- L/R States @ interface -------------------------------------------------------//
void Flux_States  (int imax, vector<vector<double> > V_Boundary, 
                   vector<vector<vector<double> > > U_cell_center, vector<vector<double> > &U_States);
//-----------------------------------------------------------------------------------------------------------------------//

//---------------------------------------- Flux -> Van Leer (FVS) -------------------------------------------------------//
void Flux_VL  (int imax, vector<vector<double> > V_Boundary, 
               vector<vector<vector<double> > > U_cell_center, vector<vector<double> > &F);
//-----------------------------------------------------------------------------------------------------------------------//

//---------------------------------------- Flux -> Roe (FDS) ------------------------------------------------------------//
void Flux_VL  (int imax, vector<vector<double> > V_Boundary, 
               vector<vector<vector<double> > > U_cell_center, vector<vector<double> > &F);
//-----------------------------------------------------------------------------------------------------------------------//

#endif