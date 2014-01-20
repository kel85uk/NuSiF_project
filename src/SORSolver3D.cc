#include "SORSolver3D.hh"

SORSolver3D::SORSolver3D ( unsigned int imax_, unsigned int jmax_, unsigned int kmax_, real omg_, real eps_, unsigned int itermax_):
                     imax(imax_), jmax(jmax_), kmax(kmax_), itermax(itermax_), checkfrequency(1), omg(omg_), eps(eps_)
{
   std::cout<<"\nCreating SOR solver\n";
}

SORSolver3D::SORSolver3D ( const FileReader & configuration ):
                     imax(configuration.getIntParameter("imax")), jmax(configuration.getIntParameter("jmax")),
                     kmax(configuration.getIntParameter("kmax")), itermax(configuration.getIntParameter("itermax")),
                     omg(configuration.getRealParameter("omg")), eps(configuration.getRealParameter("eps"))
{
   std::cout<<"\nCreating SOR solver\n";
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

inline void SORSolver3D::SetBoundary (Array<real> & p)
{
    for (unsigned int i=0; i<=imax;i++)
        for (unsigned int k=0; k<=kmax;k++) {
            p(i,0,k) = p(i,1,k);
            p(i,jmax+1,k) = p(i,jmax,k); }
    for (unsigned int j=0; j<=jmax;j++)
        for (unsigned int k=0; k<=kmax;k++) {
            p(0,j,k) = p(1,j,k);
            p(imax+1,j,k) = p(imax,j,k); }
    for (unsigned int i=0; i<=imax;i++)
        for (unsigned int j=0; j<=jmax;j++) {
            p(i,j,0) = p(i,j,1);
            p(i,j,kmax+1) = p(i,j,kmax); }
}

bool SORSolver3D::solve( StaggeredGrid3D & grid )
{
   Array<real> & p_ = grid.p();
   Array<real> & rhs_ = grid.rhs();
   real res = 0.0,resid = 1e100, dx_2 = pow(grid.dx(),-2), dy_2 = pow(grid.dy(),-2), dz_2 = pow(grid.dz(),-2);
   real numFluid_ = grid.getNumFluid();
   unsigned int iterno = 0;
   SetBoundary(p_);
// SOR iteration
   do {
      iterno++;
      for (unsigned int i=1; i<=imax;i++)
          for (unsigned int j=1; j<=jmax;j++)
              for (unsigned int k=1; k<=kmax;k++)
                  if (grid.isFluid(i,j,k))
                     p_(i,j,k) = (1.0-omg)* p_(i,j,k)
                                + omg * (  (grid.p(i,j,k,EAST) + grid.p(i,j,k,WEST)) *dx_2        
                                          +(grid.p(i,j,k,NORTH) + grid.p(i,j,k,SOUTH)) *dy_2
                                          +(grid.p(i,j,k,UP) + grid.p(i,j,k,DOWN)) *dz_2
                                          - rhs_(i,j,k))
                                      / (2.0*(dx_2 + dy_2 + dz_2));
    
      SetBoundary(p_);
      if (iterno%checkfrequency == 0 ) {
         resid =0.0;
         for (unsigned int i=1; i<=imax;i++)
	       for (unsigned int j=1; j<=jmax;j++)
	           for (unsigned int k=1; k<=kmax;k++)
	               if (grid.isFluid(i,j,k)) {
	                  res  =   (grid.p(i,j,k,EAST) + grid.p(i,j,k,WEST) - 2.0*grid.p(i,j,k,CENTER)) *dx_2
	                         + (grid.p(i,j,k,NORTH) + grid.p(i,j,k,SOUTH) - 2.0*grid.p(i,j,k,CENTER)) *dy_2
	                         + (grid.p(i,j,k,UP) + grid.p(i,j,k,DOWN) - 2.0*grid.p(i,j,k,CENTER)) *dz_2
	                         - rhs_(i,j,k) ;       
	                  resid += res*res;  }
         resid = sqrt(resid/ numFluid_);
         std::cout<<"Iteration no = "<<iterno<< "\tResidual = "<< resid<<"\n"; } }
  while (resid>=eps && iterno<itermax);   

  if (iterno<itermax) {
      std::cout<<"Solution converged in "<<iterno<< " iterations!\n";
      return true; }
  else
      std::cout<<"\nSolution not converged!\n";
  return false;
}
