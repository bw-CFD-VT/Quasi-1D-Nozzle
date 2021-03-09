// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#include "WriteFile.hpp"

using namespace std;


void norm_file(string filename, int counter, int imax, vector<double>v)
{   
    ofstream file(filename,std::ios_base::app);
    if (file.is_open())
    {
        file<<counter<<",";
        //Write Variable Outputs to Corresponding Column/Variable Label
        for (int i = 0; i<3; i++)
        {   
            file << fixed << setprecision(14) << v[i]<<",";
        }
            file << "\n";
            file.close();
    
    }
}

void mach_file(string filename, int counter, int imax, vector<double> v)
{   
    ofstream file(filename,std::ios_base::app);
    if (file.is_open())
    {
        file<<counter<<",";
        //Write Variable Outputs to Corresponding Column/Variable Label
        for (int i = 0; i<imax; i++)
        {   
            file << fixed << setprecision(14) << v[i]<<",";
        }
            file << "\n";
            file.close();
    
    }
}

void rho_file(string filename, int counter, int imax, vector<vector<vector<double> > >v)
{   
    ofstream file(filename,std::ios_base::app);
    if (file.is_open())
    {
        file<<counter<<",";
        //Write Variable Outputs to Corresponding Column/Variable Label
        for (int i = 0; i<imax; i++)
        {   
            file << fixed << setprecision(14) << v[0][i][0]<<", ";
        }
            file << "\n";
            file.close();
    
    }
}
void u_file(string filename, int counter, int imax, vector<vector<vector<double> > >v)
{   
    ofstream file(filename,std::ios_base::app);
    if (file.is_open())
    {
        file<<counter<<",";
        //Write Variable Outputs to Corresponding Column/Variable Label
        for (int i = 0; i<imax; i++)
        {   
            file << fixed << setprecision(14) << v[0][i][1]<<", ";
        }
            file << "\n";
            file.close();
    
    }
}
void press_file(string filename, int counter, int imax, vector<vector<vector<double> > >v)
{   
    ofstream file(filename,std::ios_base::app);
    if (file.is_open())
    {
        file<<counter<<",";
        //Write Variable Outputs to Corresponding Column/Variable Label
        for (int i = 0; i<imax; i++)
        {   
            file << fixed << setprecision(14) << v[0][i][2]<<", ";
        }
            file << "\n";
            file.close();
    
    }
}

void exact_file(string filename, const vector<vector<double> > v)
{   
    ofstream file(filename);
    if (file.is_open())
    {
        //Write Variable Outputs to Corresponding Column/Variable Label
        for (int i = 0; i<v[0].size(); i++)
        {   
            for (int j = 0; j<v.size(); j++)
            {
                file << fixed << setprecision(14) << v[j][i]<<",";
            }
            file << "\n";
        }
        file.close();
    }
}


void clear_existing_file()
{

    if( remove( "norm.txt" ) != 0 || remove( "mach.txt" ) != 0 || remove( "rho.txt" ) != 0 
        || remove( "u.txt" ) != 0 || remove( "press.txt" ) != 0)

    perror( "Error deleting file" );

  else

    puts( "File successfully deleted" );
}