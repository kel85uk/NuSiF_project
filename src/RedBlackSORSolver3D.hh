#ifndef RedBlackSOR_SOLVER_HH
#define RedBlackSOR_SOLVER_HH

#include "StaggeredGrid3D.hh"
#include "Solver3D.hh"
#include <omp.h>

class RedBlackSORSolver3D: public Solver3D
{
public:
   // Constructor to manually create RedBlackSORSolver3D
   RedBlackSORSolver3D ( unsigned int imax_, unsigned int jmax_, unsigned int kmax_, real omg_, real eps_, unsigned int itermax_);

   // Constructor to create a RedBlackSORSolver3D from a parsed configuration file
   RedBlackSORSolver3D ( const FileReader & configuration );


   // Solve the pressure equation on the staggered Grid3D. Returns true if converged
   bool solve( StaggeredGrid3D & Grid3D );

private:
   unsigned int imax, jmax, kmax, itermax, checkfrequency;
   real omg, eps;
   
   // copy inner points to boundary
   inline void SetBoundary (Array<real> & p);
};
#endif //RedBlackSOR_SOLVER_HH
