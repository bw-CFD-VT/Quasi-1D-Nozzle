// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#ifndef ISENTROPIC_FLOW_HPP
#define ISENTROPIC_FLOW_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window

using namespace std;

//---------------------------------------- Isentropic Flow --------------------------------------------------------------//
void Isentropic_Flow (double Mach, double &rho, double &u, double &p, double &T);
//-----------------------------------------------------------------------------------------------------------------------//

#endif