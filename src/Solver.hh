#ifndef SOLVER_HH
#define SOLVER_HH

#include "StaggeredGrid.hh"

class Solver
{
public:

   // Constructor to create a Solver
   Solver () {}

   // Solve the pressure equation on the staggered grid. Returns true if converged
   bool solve( StaggeredGrid & grid ) { return false;}

};
#endif //SOLVER_HH
