#include <iostream>
#include "FluidSimulator3D.hh"

int main()
{
   FileReader myReader;
   bool res = myReader.readFile ( "DerivativeTest3DInput.par" );
   CHECK_MSG(res, "Could not open file 'DerivativeTestInput.par' which has to be in the current directory.");
   
   real resid1 =0.0, resid2 =0.0, resid3 =0.0;   

   FluidSimulator3D simulator1(myReader);
   StaggeredGrid3D &grid1 = simulator1.grid();
   Array<real> & f1 = grid1.f();
   Array<real> & g1 = grid1.g();
   Array<real> & h1 = grid1.h();
   Array<real> & u1 = grid1.u();
   Array<real> & v1 = grid1.v();
   Array<real> & w1 = grid1.w();
   for (int i=0; i<u1.getSize(0);i++)
       for (int j=0; j<u1.getSize(1);j++) 
           for (int k=0; k<u1.getSize(1);k++) 
               u1(i,j,k) = (i+j+k)*(i+j+k);
   for (int i=0; i<v1.getSize(0);i++)
       for (int j=0; j<v1.getSize(1);j++) 
           for (int k=0; k<v1.getSize(1);k++) 
               v1(i,j,k) = i*i + j*j + k*k; 

   for (int i=0; i<w1.getSize(0);i++)
       for (int j=0; j<w1.getSize(1);j++) 
           for (int k=0; k<w1.getSize(1);k++) 
               w1(i,j,k) = i*j + j*k + k*i; 
           
   simulator1.simulateTimeStepCount();

   FluidSimulator3D simulator2(myReader);
   StaggeredGrid3D &grid2 = simulator2.grid();
   Array<real> & f2 = grid2.f();
   Array<real> & g2 = grid2.g();
   Array<real> & h2 = grid2.h();
   Array<real> & u2 = grid2.u();
   Array<real> & v2 = grid2.v();
   Array<real> & w2 = grid2.w();
   for (int i=0; i<u2.getSize(0);i++)
       for (int j=0; j<u2.getSize(1);j++) 
           for (int k=0; k<u2.getSize(1);k++) 
               w2(i,j,k) = (i+j+k)*(i+j+k);
   for (int i=0; i<v2.getSize(0);i++)
       for (int j=0; j<v2.getSize(1);j++) 
           for (int k=0; k<v2.getSize(1);k++) 
               u2(i,j,k) = i*i + j*j + k*k;

   for (int i=0; i<w2.getSize(0);i++)
       for (int j=0; j<w2.getSize(1);j++) 
           for (int k=0; k<w2.getSize(1);k++) 
               v2(i,j,k) = i*j + j*k + k*i; 
           
   simulator2.simulateTimeStepCount();

   for (int i=1; i<f1.getSize(0);i++)
       for (int j=1; j<f1.getSize(1);j++)
           for (int k=0; k<f1.getSize(1);k++) {
               resid1 += g1(i,j,k) - f2(k,i,j);
               resid2 += h1(i,j,k) - g2(k,i,j);
               resid3 += f1(i,j,k) - h2(j,i,k); }
   resid1 /= ((f1.getSize(0)-1)*(f1.getSize(1)-1));
   resid2 /= ((f1.getSize(0)-1)*(f1.getSize(1)-1));
   resid3 /= ((f1.getSize(0)-1)*(f1.getSize(1)-1));
   
   CHECK_MSG(resid1 <= 1e-8 && resid2 <= 1e-8 && resid3 <= 1e-8, "Derivative test failed!");
   std::cout<<"\nDerivative test 3D passed!\n";
   return 0;
}
