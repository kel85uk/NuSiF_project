#include "CGSolver3D.hh"

CGSolver3D::CGSolver3D ( const FileReader & configuration ):
                     imax(configuration.getIntParameter("imax")), jmax(configuration.getIntParameter("jmax")),
                     kmax(configuration.getIntParameter("kmax")), itermax(configuration.getIntParameter("itermax")),
                     eps(configuration.getRealParameter("eps")),
			   r(imax+2,jmax+2,kmax+2), Ark(imax+2,jmax+2,kmax+2), Residual(imax+2,jmax+2,kmax+2)
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

inline void CGSolver3D::SetBoundary (Array<real> & p)
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

bool CGSolver3D::solve( StaggeredGrid3D & grid )
{
	unsigned int i,j,k,numFluid_ = grid.getNumFluid();
	real dx_2,dy_2,dz_2;
	real delx = grid.dx();
	real dely = grid.dy();	
	real delz = grid.dz();	
	Array<real> &p_ = grid.p();
	Array<real> &rhs_ = grid.rhs();
	dx_2 = 1./(delx*delx);
	dy_2 = 1./(dely*dely);
	dz_2 = 1./(delz*delz);
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
      	    for (k=1;k<=kmax;k++)
			  if (grid.isFluid(i,j,k))   
		        {
                        // res = b - Ax
				Residual(i,j,k) =   (grid.p(i,j,k,EAST) + grid.p(i,j,k,WEST) - 2.0*p_(i,j,k)) *dx_2
                                          + (grid.p(i,j,k,NORTH) + grid.p(i,j,k,SOUTH) - 2.0*p_(i,j,k)) *dy_2
                                          + (grid.p(i,j,k,UP) + grid.p(i,j,k,DOWN) - 2.0*p_(i,j,k)) *dz_2
                                          - rhs_(i,j,k) ; 
				r(i,j,k) = Residual(i,j,k); 		   // p = res
				resdot += Residual(i,j,k)*Residual(i,j,k); // resdot_0
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
     			    for (k = 1; k <= kmax; ++k)
				  if (grid.isFluid(i,j,k)) {
					Ark(i,j,k) = - (r(i+1,j,k)-2*r(i,j,k)+r(i-1,j,k))*dx_2
					             - (r(i,j+1,k)-2*r(i,j,k)+r(i,j-1,k))*dy_2
					             - (r(i,j,k+1)-2*r(i,j,k)+r(i,j,k-1))*dz_2;
					rkdot += r(i,j,k)*(Ark(i,j,k));
				  }
		alphak = resdot/rkdot;
		for (i = 1; i <= imax; ++i)
			for (j = 1; j <= jmax; ++j)
     			    for (k = 1; k <= kmax; ++k)
				  if (grid.isFluid(i,j,k))
					p_(i,j,k) = p_(i,j,k) + alphak*r(i,j,k);    // Update pressure

		SetBoundary(p_);
		if (iterno%checkfrequency == 0){
		    for (i = 1; i <= imax; ++i)
			  for (j = 1; j <= jmax; ++j)
     			      for (k = 1; k <= kmax; ++k)
				    if (grid.isFluid(i,j,k))
            				Residual(i,j,k) =   (grid.p(i,j,k,EAST) + grid.p(i,j,k,WEST) - 2.0*p_(i,j,k)) *dx_2
                                                      + (grid.p(i,j,k,NORTH) + grid.p(i,j,k,SOUTH) - 2.0*p_(i,j,k)) *dy_2
                                                      + (grid.p(i,j,k,UP) + grid.p(i,j,k,DOWN) - 2.0*p_(i,j,k)) *dz_2
                                                      - rhs_(i,j,k) ; 
                                                    
		    std::cout<<"Iteration no = "<<iterno<< "\tResidual = "<< resid<<"\n";
		} 
		else
		    for (i = 1; i <= imax; ++i)
			  for (j = 1; j <= jmax; ++j)
     			      for (k = 1; k <= kmax; ++k)
				    if (grid.isFluid(i,j,k))
						Residual(i,j,k) = Residual(i,j,k) - alphak*Ark(i,j,k); // Update residual		

		resdotold = resdot;
		resdot = 0.;

		for (i = 1; i <= imax; ++i)
			for (j = 1; j <= jmax; ++j)
     			    for (k = 1; k <= kmax; ++k)
				  if (grid.isFluid(i,j,k))
					resdot += Residual(i,j,k)*Residual(i,j,k);
					
		betak = resdot/resdotold;
		
		// r_k Update search direction vector
		for (i = 1; i <= imax; ++i)
			for (j = 1; j <= jmax; ++j)
     			    for (k = 1; k <= kmax; ++k)
				  if (grid.isFluid(i,j,k))
					r(i,j,k) = Residual(i,j,k) + betak*r(i,j,k);
					
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
