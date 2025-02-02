// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#ifndef L2_NORM_HPP
#define L2_NORM_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window

using namespace std;


//------------------------- Calculate L2 Norm----------------------------------------------------------------------------//
void L2_Norm(double imax,vector<vector<double> >Residual, vector<double> &L2);
//-----------------------------------------------------------------------------------------------------------------------//


#endif
