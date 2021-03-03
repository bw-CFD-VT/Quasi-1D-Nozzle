#ifndef NORMS_HPP
#define NORMS_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window

using namespace std;


//------------------------- Calculate L Norms------------------------------------ ---------------------------------------//
void Norms (int counter,double imax, vector<vector<vector<double> > > Residual, vector<vector<vector<double> > > &L2, 
            vector<vector<vector<double> > > &L1, vector<vector<vector<double> > > &Linf);
//-----------------------------------------------------------------------------------------------------------------------//


#endif
