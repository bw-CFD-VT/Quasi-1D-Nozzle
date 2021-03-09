// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window

using namespace std;

//---------------------------------------- Boundary Conditions ----------------------------------------------------------//
void Boundary_Conditions(int imax, int Case_Flag,int ghost_cell,vector<vector<vector<double> > > V_cell_center,
                         vector<double> M_cell_center,vector<vector<double> > &V_Boundary,
                         vector<double> &M_Boundary,vector<double> &V_ghost_inflow, vector<double> &V_ghost_outflow);
//-----------------------------------------------------------------------------------------------------------------------//

#endif