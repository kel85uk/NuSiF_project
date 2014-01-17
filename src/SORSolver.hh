#ifndef SOR_SOLVER_HH
#define SOR_SOLVER_HH

#include "StaggeredGrid.hh"
#include "Solver.hh"

class SORSolver: public Solver
{
public:
   // Constructor to manually create SORSolver
   SORSolver ( unsigned int imax_, unsigned int jmax_, real omg_, real eps_, unsigned int itermax_);

   // Constructor to create a SORSolver from a parsed configuration file
   SORSolver ( const FileReader & configuration );


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
#endif //SOR_SOLVER_HH




