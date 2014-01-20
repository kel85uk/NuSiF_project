#ifndef CG_SOLVER_HH
#define CG_SOLVER_HH

#include "StaggeredGrid.hh"
#include "Solver.hh"
#include "MatrixCOO.hh"

class PCGSolver: public Solver
{
public:
   // Constructor to manually create CGSolver
   PCGSolver ( unsigned int imax_, unsigned int jmax_, real omg_, real eps_, unsigned int itermax_);

   // Constructor to create a CGSolver from a parsed configuration file
   PCGSolver ( const FileReader & configuration );


   // Solve the pressure equation on the staggered grid. Returns true if converged
   bool solve( StaggeredGrid & grid );

private:
   unsigned int imax,jmax, itermax, checkfrequency;
   real omg, eps;
   MatrixCOO Amat_;
   
   //Matrix assembly (ideally we don't need to do this every timestep!)
   void mat_assemble( StaggeredGrid& grid );
   // copy inner points to boundary
   inline void SetBoundary (Array<real> & p);   
};
#endif //CG_SOLVER_HH




