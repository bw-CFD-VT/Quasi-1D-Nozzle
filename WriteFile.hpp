#ifndef WRITE_FILE_HPP
#define WRITE_FILE_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window
#include <fstream>	//file stream library, i/o of files

using namespace std;

void norm_file(string filename, int counter, int imax, vector<vector<double> > v);
void mach_file(string filename, int counter, int imax, vector<vector<double> > v);
void rho_file(string filename, int counter, int imax, vector<vector<vector<double> > >v);
void press_file(string filename, int counter, int imax, vector<vector<vector<double> > >v);
void u_file(string filename, int counter, int imax, vector<vector<vector<double> > >v);
void exact_file(string filename, const vector<vector<double> > v);


#endif
