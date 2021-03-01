#ifndef FLUX_HPP
#define FLUX_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window

using namespace std;


//---------------------------------------- Initial Conditions -----------------------------------------------------------//
void Flux (int counter,int NI, vector<vector<vector<double>>> V_cell_center, vector<vector<vector<double>>> &V_interface);
//-----------------------------------------------------------------------------------------------------------------------//


#endif