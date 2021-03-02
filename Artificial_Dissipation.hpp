#ifndef ARTIFICIAL_DISSIPATION_HPP
#define ARTIFICIAL_DISSIPATION_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window

using namespace std;


//------------------------- Artificial Dissipation Vector, d @ iteration = counter ---------------------------------------//
void Artifical_Dissipation (double K_2, double K_4,int counter,int imax, int NI, vector<vector<double>> lambda_max, 
                            vector<vector<vector<double>>> V_cell_center, vector<vector<vector<double>>> U_cell_center,
                            vector<vector<double>> V_ghost_inflow, vector<vector<double>> V_ghost_outflow,
                            vector<vector<vector<double>>> &d);
//-----------------------------------------------------------------------------------------------------------------------//


#endif
