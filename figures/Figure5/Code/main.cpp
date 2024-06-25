
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <vector>
#include "Simulation.h"
using namespace std;

//use standard namespace
//using namespace std;



int main() 
{


       Simulation sim;
       sim.load();  
       sim.init();
       int NumIter;
       cout<< "Enter Number of Iterations:";
       cin>> NumIter;
       cout<< "NumIter "<< NumIter;
       for (int ii=1; ii<=NumIter;ii++){
      		sim.init_ic(ii);         
       		sim.run();   
     		sim.close(); 
       }
}




