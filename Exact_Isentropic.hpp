#ifndef EXACT_ISENTROPIC_HPP
#define EXACT_ISENTROPIC_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window

using namespace std;

void Isentropic_Nozzle_Exact (int imax, vector<double> x, vector<double> Nozzle_Area);

void Mach(double M_initial, double A_bar, double es, double &M);


#endif