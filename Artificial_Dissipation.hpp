#ifndef ARTIFICIAL_DISSIPATION_HPP
#define ARTIFICIAL_DISSIPATION_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window
#include <algorithm> //needed for identifying certain "characteristics" about range of elements
using namespace std;


//------------------------- Artificial Dissipation Vector, d @ iteration = counter ---------------------------------------//
void Artifical_Dissipation (double K_2, double K_4,int counter,int imax, int NI, int ghost_cell, vector<vector<double> > lambda_max, 
                            vector<vector<vector<double> > > V_cell_center, vector<vector<vector<double> > > U_cell_center,
                            vector<vector<vector<double> > > V_ghost_inflow, vector<vector<vector<double> > >U_ghost_inflow,
                            vector<vector<vector<double> > > V_ghost_outflow, vector<vector<vector<double> > >U_ghost_outflow,
                            vector<vector<vector<double> > > &d);
//-----------------------------------------------------------------------------------------------------------------------//

//------------------------- max value of pressure sensor to be used at ith cell -----------------------------------------//
double max_v (double v_im1, double v_i, double v_ip1, double v_ip2);
//-----------------------------------------------------------------------------------------------------------------------//


#endif
