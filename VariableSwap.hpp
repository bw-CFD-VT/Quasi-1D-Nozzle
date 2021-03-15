// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#ifndef VARIABLESWAP_HPP
#define VARIABLESWAP_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window

using namespace std;

//---------------------------------------- V -> U -----------------------------------------------------------------------//
void primitive_to_conserved (vector<double> V, vector<double> &U);
//-----------------------------------------------------------------------------------------------------------------------//

//---------------------------------------- U -> V -----------------------------------------------------------------------//
void conserved_to_primitive (vector<double> U, vector<double> &V);
//-----------------------------------------------------------------------------------------------------------------------//


#endif