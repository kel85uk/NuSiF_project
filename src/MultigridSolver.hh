#ifndef Multigrid_SOLVER_HH
#define Multigrid_SOLVER_HH

#include "StaggeredGrid.hh"
#include "Solver.hh"

class MultigridSolver: public Solver
{
public:
   // Constructor to manually create MultigridSolver
   MultigridSolver ( unsigned int imax_, unsigned int jmax_, real omg_, real eps_, unsigned int itermax_);

   // Constructor to create a MultigridSolver from a parsed configuration file
   MultigridSolver ( const FileReader & configuration );


   // Solve the pressure equation on the staggered grid. Returns true if converged
   bool solve( StaggeredGrid & grid );

private:
   unsigned int imax, jmax, itermax, checkfrequency, levels, itercoarse, iterpre, iterpost;
   real omg, eps;
   
   // copy inner points to boundary
   void SetBoundary (Array<real> & p);

   // Evaluate residual
   real Residual ( Array<real> & sol, Array<real> & rhs, real dx, real dy );
   
   void VCycle(unsigned int lev, Array<real> ** sol, Array<real> ** rhs, real dx, real dy);
   
   void Restrict( Array<real>& fine, Array<real>& coarse );
   void SORsweep( Array<real> & sol, Array<real> & rhs, real dx, real dy );
   void RestrictResisdual(Array<real>& sol, Array<real>& rhs, Array<real>& rhs_coarse, real dx, real dy);
   void interpolateCorrect(Array<real>& sol_fine, Array<real>& sol_coarse);
   void interpolate(Array<real>& sol_fine, Array<real>& sol_coarse);
   
   
   real rhsSum( Array<real> & rhs );

};
#endif //Multigrid_SOLVER_HH




