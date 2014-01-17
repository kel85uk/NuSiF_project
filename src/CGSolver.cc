#include "CGSolver.hh"

CGSolver::CGSolver ( unsigned int imax_, unsigned int jmax_, real omg_, real eps_, unsigned int itermax_):
                     imax(imax_), jmax(jmax_), itermax(itermax_), checkfrequency(1), omg(omg_), eps(eps_)
{
   std::cout<<"\nCreating CG solver\n";
}

CGSolver::CGSolver ( const FileReader & configuration ):
                     imax(configuration.getIntParameter("imax")), jmax(configuration.getIntParameter("jmax")),
                     itermax(configuration.getIntParameter("itermax")),eps(configuration.getRealParameter("eps"))
{
   std::cout<<"\nCreating CG solver\n";
//   CHECK_MSG((omg <= 1.9 && omg >= 1.7), "omg must be between 1.7 and 1.9");
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

// Calculate residual
inline real CGSolver::Residual ( StaggeredGrid & grid)
{
   Array<real> & rhs = grid.rhs();
   real resid=0.0, dx_2=1.0/(grid.dx()*grid.dx()), dy_2=1.0/(grid.dy()*grid.dy()), res;
   for (unsigned int i=1; i<=imax;i++)
       for (unsigned int j=1; j<=jmax;j++)
           if (grid.isFluid(i,j)) {
              res  =   (grid.p(i,j,EAST) + grid.p(i,j,WEST) - 2.0*grid.p(i,j,CENTER)) *dx_2
                     + (grid.p(i,j,NORTH) + grid.p(i,j,SOUTH) - 2.0*grid.p(i,j,CENTER)) *dy_2
                     - rhs(i,j) ;       
           resid += res*res;  }
   return (sqrt(resid/grid.getNumFluid()));
}



bool CGSolver::solve( StaggeredGrid & grid )
{
	int i,j,nfluid;
	real rdx2,rdy2;
	real add,beta_2,tol, tol0;
	real delx = grid.dx();
	real dely = grid.dy();	
	Array<real> &U = grid.p();
	Array<real> &RHS = grid.rhs();
	rdx2 = 1./delx/delx;
	rdy2 = 1./dely/dely;
	nfluid = grid.getNumFluid();
	real resid = 1e100;
	unsigned int iterno = 0;
	/** Set variables for CG */
	Array<real> Ap = U;
	Array<real> Res = U;
	Res.fill(0);
	Ap.fill(0);
	Array<real> p = U;
	p.fill(0);
	real resdot, resdotold, resdot0, pkdot, betak, alphak = 0.;
	/** Initialize CG (Calculate residual of initial condition) */
	SetBoundary(U);
	for (i=1;i<=imax;i++)
		for (j=1;j<=jmax;j++)
			if (grid.isFluid(i,j))   
		    {
				Res(i,j) = (U(i+1,j)-2*U(i,j)+U(i-1,j))*rdx2 + (U(i,j+1)-2*U(i,j)+U(i,j-1))*rdy2 - RHS(i,j); // res = b - Ax
				p(i,j) = Res(i,j); // p = res
				resdot += Res(i,j)*Res(i,j); // resdot_0
			}
	resdot0 = resdot;	
// CG iteration
   do {
      iterno++;
		/** Start CG iterations */
		pkdot = 0.;
		SetBoundary(p);
		for (i = 1; i <= imax; ++i) /* p_k Update search direction vector */
			for (j = 1; j <= jmax; ++j)
				if (grid.isFluid(i,j)){
					Ap(i,j) = -(p(i+1,j)-2*p(i,j)+p(i-1,j))*rdx2 - (p(i,j+1)-2*p(i,j)+p(i,j-1))*rdy2;
					pkdot += p(i,j)*(Ap(i,j));
				}
		alphak = resdot/pkdot;
		for (i = 1; i <= imax; ++i)
			for (j = 1; j <= jmax; ++j)
				if (grid.isFluid(i,j)){
					U(i,j) = U(i,j) + alphak*p(i,j); // Update iterate
				}
		SetBoundary(U);
		if (iterno%checkfrequency == 0){
			for (i = 1; i <= imax; ++i)
				for (j = 1; j <= jmax; ++j)
					if (grid.isFluid(i,j)){
						Res(i,j) = (U(i+1,j)-2*U(i,j)+U(i-1,j))*rdx2 + (U(i,j+1)-2*U(i,j)+U(i,j-1))*rdy2 - RHS(i,j);
					}
			std::cout<<"Iteration no = "<<iterno<< "\tResidual = "<< resid<<"\n";
		} 
		else{
			for (i = 1; i <= imax; ++i)
				for (j = 1; j <= jmax; ++j)
					if (grid.isFluid(i,j))
						Res(i,j) = Res(i,j) - alphak*Ap(i,j); // Update residual		
		}
		resdotold = resdot; // resdot_k-2
		resdot = 0.;
		for (i = 1; i <= imax; ++i) /* resdot_k-1 */
			for (j = 1; j <= jmax; ++j)
				if (grid.isFluid(i,j))
					resdot += Res(i,j)*Res(i,j);
					
		betak = resdot/resdotold; // betak = resdot_k/resdot_k-1
		
		for (i = 1; i <= imax; ++i) /* p_k Update search direction vector */
			for (j = 1; j <= jmax; ++j)
				if (grid.isFluid(i,j))
					p(i,j) = Res(i,j) + betak*p(i,j);
		resid = sqrt(resdot/nfluid);
	}
  while (resid>=eps && iterno<itermax);   

  if (iterno<itermax) {
      std::cout<<"Solution converged in "<<iterno<< " iterations!\n";
      return true; }
  else
      std::cout<<"\nSolution not converged!\n";
  return false;
}
