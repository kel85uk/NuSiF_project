#include <iostream>
#include "FluidSimulator.hh"

int main()
{
   FileReader myReader;
   bool res = myReader.readFile ( "DerivativeTestInput.par" );
   CHECK_MSG(res, "Could not open file 'DerivativeTestInput.par' which has to be in the current directory.");
   
   real resid1 =0.0, resid2 =0.0;   

   FluidSimulator simulator1(myReader);
   StaggeredGrid &grid1 = simulator1.grid();
   Array<real> & f1 = grid1.f();
   Array<real> & g1 = grid1.g();
   Array<real> & u1 = grid1.u();
   Array<real> & v1 = grid1.v();
   for (int i=0; i<u1.getSize(0);i++)
       for (int j=0; j<u1.getSize(1);j++) 
           u1(i,j) = (i+j)*(i+j);
   for (int i=0; i<v1.getSize(0);i++)
       for (int j=0; j<v1.getSize(1);j++) 
           v1(i,j) = i*i + j*j; 
           
   simulator1.simulateTimeStepCount();

   FluidSimulator simulator2(myReader);
   StaggeredGrid &grid2 = simulator2.grid();
   Array<real> & f2 = grid2.f();
   Array<real> & g2 = grid2.g();
   Array<real> & u2 = grid2.u();
   Array<real> & v2 = grid2.v();
   for (int i=0; i<u2.getSize(0);i++)
       for (int j=0; j<u2.getSize(1);j++) 
           u2(i,j) = i*i + j*j;
   for (int i=0; i<v2.getSize(0);i++)
       for (int j=0; j<v2.getSize(1);j++) 
           v2(i,j) = (i+j)*(i+j); 
   simulator2.simulateTimeStepCount();

   for (int i=1; i<f1.getSize(0);i++)
       for (int j=1; j<f1.getSize(1);j++) {
           resid1 += f1(i,j) - g2(j,i);
           resid2 += g1(i,j) - f2(j,i); }
           
   resid1 /= ((f1.getSize(0)-1)*(f1.getSize(1)-1));
   resid2 /= ((f1.getSize(0)-1)*(f1.getSize(1)-1));
   
   CHECK_MSG(resid1 <= 1e-8 && resid2 <= 1e-8, "Derivative test failed!");
   std::cout<<"\nDerivative test passed!\n";
   return 0;
}
