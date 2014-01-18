#ifndef CG_Solver3D_HH
#define CG_Solver3D_HH

#include "StaggeredGrid3D.hh"
#include "Solver3D.hh"

class CGSolver3D: public Solver3D
{
public:

   // Constructor to create a CGSolver3D from a parsed configuration file
   CGSolver3D ( const FileReader & configuration );

   // Solve the pressure equation on the staggered grid. Returns true if converged
   bool solve( StaggeredGrid3D & grid );

private:
   unsigned int imax, jmax, kmax, itermax, checkfrequency;
   real eps;
   Array<real> r;
   Array<real> Ark;
   Array<real> Residual;
   
   // copy inner points to boundary
   inline void SetBoundary (Array<real> & p);

};
#endif //CG_Solver3D_HH
