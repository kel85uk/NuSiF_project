#include <iostream>
#include <cmath>
#include "PCGSolver.hh"

// 2D Setup 1:
void initGridSetup1( StaggeredGrid & grid )
{
   Array<real> & rhs = grid.rhs();
   Array<real> & p = grid.p();
   for (int i=0; i<rhs.getSize(0);i++)
       for (int j=0; j<rhs.getSize(1);j++) {
   //    - grid.p   : init with random values
           p(i,j) =  rand() % 10000-20000;
   //    - grid.rhs : init with zero
           rhs(i,j) = 0; }
}

// 2D Setup 2:
void initGridSetup2( StaggeredGrid & grid )
{
   Array<real> & rhs = grid.rhs();
   Array<real> & p = grid.p();
   real rhs_sum = 0.0;
   for (int i=1; i<rhs.getSize(0)-1;i++)
       for (int j=1; j<rhs.getSize(1)-1;j++) {
   //    - grid.p   : init with random values
           p(i,j) =  rand() % 10000-20000;
   //    - grid.rhs : f(x,y) = sin(2 * x * \pi)
           rhs(i,j) = sin(2.0*M_PI*(i-0.5)*grid.dx());
           rhs_sum += rhs(i,j); }

   rhs_sum /= ((rhs.getSize(0)-2)*(rhs.getSize(1)-2));

   for (int i=1; i<rhs.getSize(0)-1;i++)
       for (int j=1; j<rhs.getSize(1)-1;j++)
           rhs(i,j) -= rhs_sum;

}

int main()
{
   bool converged;
   std::cout<<"\nReading file ";
   FileReader configuration;
   configuration.readFile("poisson.par");

   std::cout<<"\nMaking Staggered Grid\n";
   StaggeredGrid grid(configuration);

   std::cout<<"\nCreating Solver\n";
   PCGSolver solver(configuration);

   std::cout<<"\nInitializing Grid Setup 1\n";
   initGridSetup1( grid );
   std::cout<<"\nRunning Solver for Grid Setup 1\n";
   converged = solver.solve(grid);
   CHECK_MSG(converged , "Grid Setup 1 not converged");
//   grid.p().print();
   std::cout<<"\nInitializing Grid Setup 2\n";
   initGridSetup2( grid );
   std::cout<<"\nRunning Solver for Grid Setup 2\n";
   converged = solver.solve(grid);
   CHECK_MSG(converged, "Grid Setup 2 not converged");
//   grid.p().print();
   std::cout<<"\nCG Test 2D passed\n";
   return 0;
}
