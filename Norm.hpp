#ifndef NORM_HPP
#define NORM_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window

using namespace std;


//------------------------- Calculate L2 Norm----------------------------------------------------------------------------//
void Norm(int counter,double imax,vector<vector<vector<double> > > Residual, vector<vector<double> > &L2);
//-----------------------------------------------------------------------------------------------------------------------//


#endif
