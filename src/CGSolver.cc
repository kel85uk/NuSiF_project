#include "CGSolver.hh"

CGSolver::CGSolver ( const FileReader & configuration ):
                     imax(configuration.getIntParameter("imax")), jmax(configuration.getIntParameter("jmax")),
                     itermax(configuration.getIntParameter("itermax")),eps(configuration.getRealParameter("eps")),
			   r(imax+2,jmax+2), Ark(imax+2,jmax+2), Residual(imax+2,jmax+2)
{
   std::cout<<"\nCreating CG solver\n";
   CHECK_MSG((eps > 0), "eps must be greater than 0");
   if (!configuration.find ("checkfrequency")) {
       WARN("\nCheck Frequency not specified\n");
      checkfrequency = 1; }
   else  
      checkfrequency = configuration.getIntParameter("checkfrequency");
   CHECK_MSG((checkfrequency > 0), "check frequency must be greater than 0");
   CHECK_MSG((itermax > 0), "Max number of iterations (itermax) must be greater than 0");
}

inline void CGSolver::SetBoundary ( Array<real> & p )
{
   for (unsigned int i=0; i<=imax;i++) {
        p(i,0) = p(i,1);
        p(i,jmax + 1) = p(i,jmax); }
   for (unsigned int j=0; j<=jmax;j++) {
        p(0,j) = p(1,j);
        p(imax + 1,j) = p(imax,j); }
}

bool CGSolver::solve( StaggeredGrid & grid )
{
	unsigned int i,j,numFluid_ = grid.getNumFluid();
	real dx_2,dy_2;
	real delx = grid.dx();
	real dely = grid.dy();	
	Array<real> &p_ = grid.p();
	Array<real> &rhs_ = grid.rhs();
	dx_2 = 1./(delx*delx);
	dy_2 = 1./(dely*dely);
	real resid = 1e100;
	unsigned int iterno = 0;

	/** Set variables for CG */
	r.fill(0);
	Ark.fill(0);
	Residual.fill(0);
	real resdot = 0., resdotold, rkdot, betak, alphak = 0.;

	/** Initialize CG (Calculate residual of initial condition) */
	SetBoundary(p_);
	for (i=1;i<=imax;i++)
		for (j=1;j<=jmax;j++)
			if (grid.isFluid(i,j))   
		    {
                        // res = b - Ax
				Residual(i,j) =   (grid.p(i,j,EAST) + grid.p(i,j,WEST) - 2.0*p_(i,j)) *dx_2
                                        + (grid.p(i,j,NORTH) + grid.p(i,j,SOUTH) - 2.0*p_(i,j)) *dy_2
                                        - rhs_(i,j) ; 
				r(i,j) = Residual(i,j); 		   // p = res
				resdot += Residual(i,j)*Residual(i,j); // resdot_0
			}
// CG iteration
   do {
            iterno++;
		// Start CG iterations
		rkdot = 0.0;

		SetBoundary(r);
		
		// r_k Update search direction vector
		for (i = 1; i <= imax; ++i)
			for (j = 1; j <= jmax; ++j)
				if (grid.isFluid(i,j)){
					Ark(i,j) = -(r(i+1,j)-2*r(i,j)+r(i-1,j))*dx_2 - (r(i,j+1)-2*r(i,j)+r(i,j-1))*dy_2;
					rkdot += r(i,j)*(Ark(i,j));
				}
		alphak = resdot/rkdot;
		for (i = 1; i <= imax; ++i)
			for (j = 1; j <= jmax; ++j)
				if (grid.isFluid(i,j)){
					p_(i,j) = p_(i,j) + alphak*r(i,j);    // Update pressure
				}
		SetBoundary(p_);
		if (iterno%checkfrequency == 0){
			for (i = 1; i <= imax; ++i)
				for (j = 1; j <= jmax; ++j)
					if (grid.isFluid(i,j))
            				Residual(i,j) =   (grid.p(i,j,EAST) + grid.p(i,j,WEST) - 2.0*p_(i,j)) *dx_2
                                                    + (grid.p(i,j,NORTH) + grid.p(i,j,SOUTH) - 2.0*p_(i,j)) *dy_2
                                                    - rhs_(i,j) ; 
                                                    
			std::cout<<"Iteration no = "<<iterno<< "\tResidual = "<< resid<<"\n";
		} 
		else{
			for (i = 1; i <= imax; ++i)
				for (j = 1; j <= jmax; ++j)
					if (grid.isFluid(i,j))
						Residual(i,j) = Residual(i,j) - alphak*Ark(i,j); // Update residual		
		}

		resdotold = resdot;
		resdot = 0.;

		for (i = 1; i <= imax; ++i)
			for (j = 1; j <= jmax; ++j)
				if (grid.isFluid(i,j))
					resdot += Residual(i,j)*Residual(i,j);
					
		betak = resdot/resdotold;
		
		// r_k Update search direction vector
		for (i = 1; i <= imax; ++i)
			for (j = 1; j <= jmax; ++j)
				if (grid.isFluid(i,j))
					r(i,j) = Residual(i,j) + betak*r(i,j);
					
		resid = sqrt(resdot/numFluid_);
	}
  while (resid>=eps && iterno<itermax);   

  if (iterno<itermax) {
      std::cout<<"Solution converged in "<<iterno<< " iterations!\n";
      return true; }
  else
      std::cout<<"\nSolution not converged!\n";
  return false;
}
