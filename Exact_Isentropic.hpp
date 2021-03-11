// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#ifndef EXACT_ISENTROPIC_HPP
#define EXACT_ISENTROPIC_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window

using namespace std;

//--------------------------------------- Calculate Exact Solution (Isentropic) -----------------------------------------//
void Exact_Solution(int imax, vector<double> x_cell_center, vector<double> Area_cell_center,
                    vector<double> M_exact, vector<double> rho_exact, vector<double> u_vel_exact,
                    vector<double> p_exact, vector<double> T_exact, vector<vector<double> > V_exact, vector<vector<double> >&U_exact);
//-----------------------------------------------------------------------------------------------------------------------//


//------------------------------------ Newton Iteration to Calculate Exact Mach -----------------------------------------//
void Mach_Exact(double M_initial, double A_bar, double es, double &M);
//-----------------------------------------------------------------------------------------------------------------------//


#endif