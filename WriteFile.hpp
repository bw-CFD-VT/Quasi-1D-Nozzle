// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#ifndef WRITE_FILE_HPP
#define WRITE_FILE_HPP

#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>	//std io stream library, keyboard and terminal window
#include <fstream>	//file stream library, i/o of files
#include <iomanip>

using namespace std;

void residual_norm_file(string filename, int counter, int imax, vector<double>v);
void error_norm_file(string filename, int imax, vector<double> v);
void mach_file(string filename, int counter, int imax, vector<double> v);
void prim_variable_file(string filename, int variable, int counter, int imax, vector<vector<vector<double> > >v);
void exact_file(string filename, vector<vector<double> > v);
void clear_existing_file();

#endif
