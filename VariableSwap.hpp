#ifndef VARIABLESWAP_HPP
#define VARIABLESWAP_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window

using namespace std;

//---------------------------------------- V -> U -----------------------------------------------------------------------//
void primative_to_conserved (int counter, vector<vector<vector<double> > > V, vector<vector<vector<double> > > &U);
//-----------------------------------------------------------------------------------------------------------------------//

//---------------------------------------- U -> V -----------------------------------------------------------------------//
void conserved_to_primative (int counter, vector<vector<vector<double> > > U, vector<vector<vector<double> > > &V);
//-----------------------------------------------------------------------------------------------------------------------//


#endif