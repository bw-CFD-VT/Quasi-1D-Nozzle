#ifndef INITIAL_CONDITIONS_HPP
#define INITIAL_CONDITIONS_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window

using namespace std;

//---------------------------------------- Initial Conditions -----------------------------------------------------------//
void Initial_Conditions (int imax,int NI,vector<double> x_cell_center, 
                         vector<vector<vector<double> > > &V_cell_center, vector<vector<double> > &M_cell_center);
//-----------------------------------------------------------------------------------------------------------------------//


#endif