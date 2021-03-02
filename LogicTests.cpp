#include <vector>  //needed for vectors
#include <cmath>   //needed for math operations
#include <iostream>
using namespace std;

int main()
{

    vector<vector<double>> x;

    x.resize(1);

    x[0] ={1,2,3};

    cout<<x[0][0]<<", "<<x[0][1]<<", "<<x[0][2]<<endl;

    x.resize(2);
    x[1].push_back(4);
  cout<<x[0][0]<<", "<<x[0][1]<<", "<<x[0][2]<<endl;

    cout<<x[1][0]<<endl;


}