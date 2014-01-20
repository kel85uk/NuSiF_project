#ifndef CG_SOLVER_HH
#define CG_SOLVER_HH

#include "StaggeredGrid.hh"
#include "Solver.hh"

class CGSolver: public Solver
{
public:

   // Constructor to create a CGSolver from a parsed configuration file
   CGSolver ( const FileReader & configuration );

   // Solve the pressure equation on the staggered grid. Returns true if converged
   bool solve( StaggeredGrid & grid );

private:
   unsigned int imax,jmax, itermax, checkfrequency;
   real eps;
   Array<real> r;
   Array<real> Ark;
   Array<real> Residual;
   
   // copy inner points to boundary
   inline void SetBoundary (Array<real> & p);

};
#endif //CG_SOLVER_HH
