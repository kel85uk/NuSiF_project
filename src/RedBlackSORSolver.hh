#ifndef RedBlackSOR_SOLVER_HH
#define RedBlackSOR_SOLVER_HH

#include "StaggeredGrid.hh"
#include "Solver.hh"
#include <omp.h>

class RedBlackSORSolver: public Solver
{
public:
   // Constructor to manually create RedBlackSORSolver
   RedBlackSORSolver ( unsigned int imax_, unsigned int jmax_, real omg_, real eps_, unsigned int itermax_);

   // Constructor to create a RedBlackSORSolver from a parsed configuration file
   RedBlackSORSolver ( const FileReader & configuration );


   // Solve the pressure equation on the staggered grid. Returns true if converged
   bool solve( StaggeredGrid & grid );

private:
   unsigned int imax,jmax, itermax, checkfrequency;
   real omg, eps;
   
   // copy inner points to boundary
   inline void SetBoundary (Array<real> & p);

   // Evaluate residual
   inline real Residual (StaggeredGrid & grid);   
};
#endif //RedBlackSOR_SOLVER_HH




