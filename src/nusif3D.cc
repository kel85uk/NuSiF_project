#include <iostream>
#include "FluidSimulator3D.hh"

int main( int argc, char* argv[])
{

   CHECK_MSG(argc>1, "Filename missing");
   std::cout<<"\nReading file " << argv[1] <<"\n" ;
   FileReader myReader;
   myReader.readFile(argv[1]);

   FluidSimulator3D simulator(myReader);
   StaggeredGrid3D &grid = simulator.grid();
   Array<real> & rhs = grid.rhs();
   for (int i=0; i<rhs.getSize(0);i++)
       for (int j=0; j<rhs.getSize(1);j++) 
           for (int k=0; k<rhs.getSize(2);k++) 
               rhs(i,j,k) = 0;

   simulator.simulateTimeStepCount();   
   return 0;
}
