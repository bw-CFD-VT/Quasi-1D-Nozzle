#include <iostream>	//std io stream library, keyboard and terminal window
#include <fstream>	//file stream library, i/o of files
#include <string>   //needed for strings
#include <iomanip>  //fieldwidth manipulators allowed to be used

// void csv_file(string filename, vector<double> v1,vector<double> v2, vector<double> v3, vector<double> v4)
// {   
//     ofstream file(filename);
//     if (file.is_open())
//     {
//         //Write Column Labels as Variables
//         for (int i = 0; i<v1.size(); k++)
//         {
//             file << v1[k]<<",";
//         }

//         file << "\n";


//         //Write Variable Outputs to Corresponding Column/Variable Label
//         for (int i = 0; i<v[0].size(); i++)
//         {   
//             for (int j = 0; j<v.size(); j++)
//             {
//                 file << fixed << setprecision(14) << v[j][i]<<",";
//             }
//             file << "\n";
//         }
//         file.close();
//     }
// }