#include <iostream>
#include "FluidSimulator.hh"

int main( int argc, char* argv[])
{

   CHECK_MSG(argc>1, "Filename missing");
   std::cout<<"\nReading file " << argv[1] <<"\n" ;
   FileReader myReader;
   myReader.readFile(argv[1]);

   FluidSimulator simulator(myReader);
   StaggeredGrid &grid = simulator.grid();
   Array<real> & rhs = grid.rhs();
   for (int i=0; i<rhs.getSize(0);i++)
       for (int j=0; j<rhs.getSize(1);j++) 
           rhs(i,j) = 0;

   simulator.simulateTimeStepCount();   
   return 0;
}
