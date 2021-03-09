// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#ifndef SOURCE_TERM_HPP
#define SOURCE_TERM_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window

using namespace std;


//------------------------- Calculate Source Term -----------------------------------------------------------------------//
void Source_Term (int counter,double imax,double dx, vector<double> Area_interface,
                  vector<vector<vector<double> > > V_cell_center,vector<vector<double> > &SourceTerm);
//-----------------------------------------------------------------------------------------------------------------------//


#endif