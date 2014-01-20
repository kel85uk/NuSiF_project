#ifndef SOR_SOLVER_3D_HH
#define SOR_SOLVER_3D_HH

#include "StaggeredGrid3D.hh"
#include "Solver3D.hh"

class SORSolver3D: public Solver3D
{
public:
   // Constructor to manually create SORSolver
   SORSolver3D ( unsigned int imax_, unsigned int jmax_, unsigned int kmax_, real omg_, real eps_, unsigned int itermax_);

   // Constructor to create a SORSolver from a parsed configuration file
   SORSolver3D ( const FileReader & configuration );


   // Solve the pressure equation on the staggered grid. Returns true if converged
   bool solve( StaggeredGrid3D & grid );

private:
   unsigned int imax, jmax, kmax, itermax, checkfrequency;
   real omg, eps;
   
   // copy inner points to boundary
   inline void SetBoundary (Array<real> & p);

};

#endif //SOR_SOLVER_3D_HH
