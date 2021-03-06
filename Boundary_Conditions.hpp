#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window

using namespace std;

//---------------------------------------- Boundary Conditions ----------------------------------------------------------//
void Boundary_Conditions(int counter,int Case_Flag,int ghost_cell,int imax,int NI,vector<vector<vector<double> > > V_cell_center,
                         vector<vector<double> > M_cell_center,vector<vector<vector<double> > > &V_Boundary,
                         vector<vector<double> > &M_Boundary,vector<vector<vector<double> > > &V_ghost_inflow,
                         vector<vector<vector<double> > > &V_ghost_outflow);
//-----------------------------------------------------------------------------------------------------------------------//

#endif