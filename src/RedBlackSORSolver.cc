#include "RedBlackSORSolver.hh"

RedBlackSORSolver::RedBlackSORSolver ( unsigned int imax_, unsigned int jmax_, real omg_, real eps_, unsigned int itermax_):
                     imax(imax_), jmax(jmax_), itermax(itermax_), checkfrequency(1), omg(omg_), eps(eps_)
{
   std::cout<<"\nCreating RedBlackSOR solver\n";
}

RedBlackSORSolver::RedBlackSORSolver ( const FileReader & configuration ):
                     imax(configuration.getIntParameter("imax")), jmax(configuration.getIntParameter("jmax")),
                     itermax(configuration.getIntParameter("itermax")),
                     omg(configuration.getRealParameter("omg")), eps(configuration.getRealParameter("eps"))
{
   std::cout<<"\nCreating RedBlackSOR solver\n";
   CHECK_MSG((omg <= 1.9 && omg >= 1.7), "omg must be between 1.7 and 1.9");
   CHECK_MSG((eps > 0), "eps must be greater than 0");
   if (!configuration.find ("checkfrequency")) {
       WARN("\nCheck Frequency not specified\n");
      checkfrequency = 1; }
   else  
      checkfrequency = configuration.getIntParameter("checkfrequency");
   CHECK_MSG((checkfrequency > 0), "check frequency must be greater than 0");
   CHECK_MSG((itermax > 0), "Max number of iterations (itermax) must be greater than 0");
}

inline void RedBlackSORSolver::SetBoundary ( Array<real> & p )
{
   for (unsigned int i=0; i<=imax;i++) {
        p(i,0) = p(i,1);
        p(i,jmax + 1) = p(i,jmax); }
   for (unsigned int j=0; j<=jmax;j++) {
        p(0,j) = p(1,j);
        p(imax + 1,j) = p(imax,j); }
}

bool RedBlackSORSolver::solve( StaggeredGrid & grid )
{
   Array<real> & p_ = grid.p();
   Array<real> & rhs_ = grid.rhs();
   real res =0.0, resid = 1e100, dx_2 = pow(grid.dx(),-2), dy_2 = pow(grid.dy(),-2);
   real numFluid = grid.getNumFluid();
   unsigned int iterno = 0;
   unsigned int i,j;
   SetBoundary(p_);
// RedBlackSOR iteration
   do {
      iterno++;
      #pragma omp parallel for private(i)
      for (j=1; j<=jmax;j++)
          for (i=1+(j+1)%2;i<=imax;i+=2)
              if (grid.isFluid(i,j))
                 p_(i,j) = (1.0-omg)* p_(i,j)
                        + omg * (  (grid.p(i,j,EAST)+grid.p(i,j,WEST)) *dx_2
                                  +(grid.p(i,j,NORTH)+grid.p(i,j,SOUTH)) *dy_2
                                  - rhs_(i,j))
                              / (2.0*(dx_2 + dy_2));
	  #pragma omp parallel for private(i)
      for (j=1; j<=jmax;j++)
          for (i=1+(j)%2;i<=imax;i+=2)
              if (grid.isFluid(i,j))
                 p_(i,j) = (1.0-omg)* p_(i,j)
                        + omg * (  (grid.p(i,j,EAST)+grid.p(i,j,WEST)) *dx_2
                                  +(grid.p(i,j,NORTH)+grid.p(i,j,SOUTH)) *dy_2
                                  - rhs_(i,j))
                              / (2.0*(dx_2 + dy_2));
    
      SetBoundary(p_);
      
      // Calculate residual
      if (iterno%checkfrequency == 0 ) {
          res = 0.0, resid = 0.0;
          for (i=1; i<=imax;i++)
              for (j=1; j<=jmax;j++)
                  if (grid.isFluid(i,j)) {
                      res  =   (grid.p(i,j,EAST) + grid.p(i,j,WEST) - 2.0*p_(i,j)) *dx_2
                             + (grid.p(i,j,NORTH) + grid.p(i,j,SOUTH) - 2.0*p_(i,j)) *dy_2
                             - rhs_(i,j) ;       
                      resid += res*res;  }
         resid = sqrt(resid/numFluid);
         std::cout<<"Iteration no = "<<iterno<< "\tResidual = "<< resid<<"\n"; } }
  while (resid>=eps && iterno<itermax);   

  if (iterno<itermax) {
      std::cout<<"Solution converged in "<<iterno<< " iterations!\n";
      return true; }
  else
      std::cout<<"\nSolution not converged!\n";
  return false;
}
