#ifndef INITIAL_BOUNDARY_CONDITIONS_HPP
#define INITIAL_BOUNDARY_CONDITIONS_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window

using namespace std;


//---------------------------------------- Initial Conditions -----------------------------------------------------------//
void Initial_Conditions (int imax,int NI,vector<double> x_cell_center, vector<double> x_interface,
                         vector<vector<vector<double>>> &V_cell_center, vector<vector<double>> &M_cell_center);
//-----------------------------------------------------------------------------------------------------------------------//


//---------------------------------------- Boundary Conditions ----------------------------------------------------------//
void Boundary_Conditions(int counter, int ghost_cell,int imax,int NI,vector<double> x_cell_center, 
                         vector<double> x_interface,vector<vector<vector<double>>> V_cell_center,
                         vector<vector<double>> M_cell_center,vector<vector<vector<double>>> &V_interface,
                         vector<vector<double>> &M_interface);
//-----------------------------------------------------------------------------------------------------------------------//


//---------------------------------------- Isentropic Flow --------------------------------------------------------------//
void Isentropic_Relationships (double Mach, double &rho, double &u, double &p, double &T);
//-----------------------------------------------------------------------------------------------------------------------//

#endif