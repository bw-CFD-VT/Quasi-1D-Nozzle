// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#ifndef ARTIFICIAL_DISSIPATION_HPP
#define ARTIFICIAL_DISSIPATION_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window
#include <algorithm> //needed for identifying certain "characteristics" about range of elements
using namespace std;


//------------------------- Artificial Dissipation Vector ---------------------------------------------------------------//
void Artifical_Dissipation (double K_2, double K_4,int imax, int ghost_cell, vector<double>lambda_max, 
                            vector<vector<vector<double> > > V_cell_center, vector<vector<vector<double> > > U_cell_center,
                            vector<double> V_ghost_inflow, vector<double> U_ghost_inflow,
                            vector<double> V_ghost_outflow, vector<double> U_ghost_outflow,
                            vector<vector<double> >&d);
//-----------------------------------------------------------------------------------------------------------------------//

//------------------------- Pressure Sensor -> Determine max(v(i-1),v(i),v(i+1),v(i+2)) ---------------------------------//
double max_v (double v_im1, double v_i, double v_ip1, double v_ip2);
//-----------------------------------------------------------------------------------------------------------------------//


#endif
