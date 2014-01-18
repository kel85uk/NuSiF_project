#ifndef SOLVER_HH
#define SOLVER_HH

#include "StaggeredGrid3D.hh"

class Solver3D
{
public:

   // Constructor to create a Solver
   Solver3D () {}

   // Solve the pressure equation on the staggered grid. Returns true if converged
   virtual bool solve( StaggeredGrid3D & grid ) { return false;}

};
#endif //SOLVER_HH
