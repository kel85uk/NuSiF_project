#include <iostream>
#include <cmath>
#include "RedBlackSORSolver3D.hh"

// 3D Setup 1:
void initGridSetup1( StaggeredGrid3D & grid )
{
   Array<real> & rhs = grid.rhs();
   Array<real> & p = grid.p();
   for (int i=0; i<rhs.getSize(0);i++)
       for (int j=0; j<rhs.getSize(1);j++) 
           for (int k=0; k<rhs.getSize(2);k++) {
       //    - grid.p   : init with random values
               p(i,j,k) =  rand() % 100-200;
       //    - grid.rhs : init with zero
               rhs(i,j,k) = 0; }
}

// 3D Setup 2:
void initGridSetup2( StaggeredGrid3D & grid )
{
   Array<real> & rhs = grid.rhs();
   Array<real> & p = grid.p();
   real rhs_sum = 0.0;
   for (int i=1; i<rhs.getSize(0)-1;i++)
       for (int j=1; j<rhs.getSize(1)-1;j++) 
           for (int k=1; k<rhs.getSize(2)-1;k++) {
       //    - grid.p   : init with random values
               p(i,j,k) =  rand() % 100-200;
       //    - grid.rhs : f(x,y,z) = sin(2 * x * \pi)
               rhs(i,j,k) = sin(2.0*M_PI*(i-0.5)*grid.dx());
               rhs_sum += rhs(i,j,k); }

   rhs_sum /= ((rhs.getSize(0)-2)*(rhs.getSize(1)-2)*(rhs.getSize(2)-2));

   for (int i=1; i<rhs.getSize(0)-1;i++)
       for (int j=1; j<rhs.getSize(1)-1;j++)
           for (int k=1; k<rhs.getSize(2)-1;k++)
           rhs(i,j,k) -= rhs_sum;

}
int main()
{
   bool converged;
   std::cout<<"\nReading file ";

   FileReader configuration3D;
   configuration3D.readFile("poisson3D.par");

   std::cout<<"\nMaking Staggered Grid\n";
   StaggeredGrid3D grid3D(configuration3D);

   std::cout<<"\nCreating Solver\n";
   RedBlackSORSolver3D solver3D(configuration3D);

   std::cout<<"\nInitializing 3DGrid Setup 1\n";
   initGridSetup1( grid3D );
   std::cout<<"\nRunning Solver for 3D Grid Setup 1\n";
   converged = solver3D.solve(grid3D);
   CHECK_MSG(converged , "3D Grid Setup 1 not converged");

   std::cout<<"\nInitializing 3D Grid Setup 2\n";
   initGridSetup2( grid3D );
   std::cout<<"\nRunning Solver for 3D Grid Setup 2\n";
   converged = solver3D.solve(grid3D);
   CHECK_MSG(converged, "3D Grid Setup 2 not converged");
   
   std::cout<<"\nRedBlackSOR Test 3D passed\n";
   return 0;
}
