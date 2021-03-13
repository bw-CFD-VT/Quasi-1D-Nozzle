// AOE 6145
// Homework 2: Quasi-1D Nozzle FVM Code
// Brendan Walsh (PID: bwalsh4)

#include "WriteFile.hpp"
#include "Constants.hpp"


using namespace std;

void residual_norm_file(string filename, int counter, int imax, vector<double>v, string grid_ID)
{   
    ofstream file(filename+grid_ID,std::ios_base::app);
    if (file.is_open())
    {
        file<<counter<<",";
        for (int i = 0; i<3; i++)
        {   
            file << fixed << setprecision(14) << v[i]<<",";
        }
            file << "\n";
            file.close();
    
    }
}

void error_norm_file(string filename, int imax, vector<double>v, string grid_ID)
{   
    ofstream file(filename+grid_ID,std::ios_base::app);
    if (file.is_open())
    {
        for (int i = 0; i<3; i++)
        {   
            file << fixed << setprecision(14) << v[i]<<",";
        }
            file << "\n";
            file.close();
    }
}

void mach_file(string filename, int counter, int imax, vector<double> v, string grid_ID)
{   
    ofstream file(filename+grid_ID,std::ios_base::app);
    if (file.is_open())
    {
        file<<counter<<",";
        for (int i = 0; i<imax; i++)
        {   
            file << fixed << setprecision(14) << v[i]<<",";
        }
            file << "\n";
            file.close();
    
    }
}

void prim_variable_file(string filename, int variable, int counter, int imax, vector<vector<vector<double> > >v, string grid_ID)
{   
    if (variable == 3)
    {
        vector<vector<vector<double> > > T(1,vector<vector<double> >(imax,vector<double>(1,0)));
        for (int i = 0; i<imax; i++)
        {
            T[0][i][0] = v[0][i][2]/v[0][i][0]/R;
        }

        variable = 0;
        v[0] = T[0];
    }

    ofstream file(filename+grid_ID,std::ios_base::app);
    if (file.is_open())
    {
        file<<counter<<",";
        for (int i = 0; i<imax; i++)
        {   
            file << fixed << setprecision(14) << v[0][i][variable]<<", ";
        }
            file << "\n";
            file.close();
    
    }
}

void exact_file(string filename, const vector<vector<double> > v, string grid_ID)
{   
    ofstream file(filename+grid_ID);
    if (file.is_open())
    {
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

void clear_existing_file(string grid_ID)
{
  string Exact_File = exactfile + grid_ID;
  string Residual_File = residualnormfile + grid_ID;
  string DE_U_File = errornormfile_U + grid_ID;
  string DE_V_File = errornormfile_V + grid_ID;
  string Mach_file = machfile + grid_ID;
  string rho_file = rhofile + grid_ID;
  string u_file = ufile + grid_ID;
  string press_file = pressfile + grid_ID;
  string temp_file = tempfile + grid_ID;

    if( remove( Exact_File.c_str() )    != 0) perror( "Error deleting file: 1" );
    if( remove( Residual_File.c_str() ) != 0) perror( "Error deleting file: 2" );
    if( remove( DE_U_File.c_str() )     != 0) perror( "Error deleting file: 3" );
    if( remove( DE_V_File.c_str() )     != 0) perror( "Error deleting file: 4" );
    if( remove( Mach_file.c_str() )     != 0) perror( "Error deleting file: 5" );
    if( remove( rho_file.c_str() )      != 0) perror( "Error deleting file: 6" );
    if( remove( u_file.c_str() )        != 0) perror( "Error deleting file: 7" );
    if( remove( temp_file.c_str() )     !=0)  perror( "Error deleting file: 8" );


    
}