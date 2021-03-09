// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#ifndef EXACT_ISENTROPIC_HPP
#define EXACT_ISENTROPIC_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window

using namespace std;

void Isentropic_Nozzle_Exact (double imax, vector<double> x, vector<double> Nozzle_Area, 
                              vector<double> &M_exact, vector<double> &rho_exact, vector<double> &u_exact,
                              vector<double> &p_exact, vector<double> &T_exact);

void Mach(double M_initial, double A_bar, double es, double &M);


#endif