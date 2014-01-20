#include "PCGSolver.hh"

PCGSolver::PCGSolver ( unsigned int imax_, unsigned int jmax_, real omg_, real eps_, unsigned int itermax_):
                     imax(imax_), jmax(jmax_), itermax(itermax_), checkfrequency(1), omg(omg_), eps(eps_)
{
   std::cout<<"\nCreating CG solver\n";
}

PCGSolver::PCGSolver ( const FileReader & configuration ):
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

inline void PCGSolver::SetBoundary ( Array<real> & p )
{
   for (unsigned int i=0; i<=imax;i++) {
        p(i,0) = p(i,1);
        p(i,jmax + 1) = p(i,jmax); }
   for (unsigned int j=0; j<=jmax;j++) {
        p(0,j) = p(1,j);
        p(imax + 1,j) = p(imax,j); }
}

bool PCGSolver::solve( StaggeredGrid & grid )
{
	std::cout << "assembling matrix" << std::endl;
	PCGSolver::mat_assemble(grid);
	unsigned int iterno = 0;
	unsigned int nfluid = grid.getNumFluid();
	int NxG = grid.p().getSize(0);
	int NyG = grid.p().getSize(1);
	Array<real> SOL = grid.p();
	Array<real> b = grid.rhs();
	real resid = 1e100;
/** Set variables for CG */
	Array<real> q = SOL;
	Array<real> Res = SOL;
	Res.fill(0);
	q.fill(0);
	Array<real> d = SOL;
	d.fill(0);
	real p0 = sqrt(SOL.dotNC(SOL));
	real dnew, dold, d0, betak, alphak = 0.;
	std::cout << "Start PCG loop" << std::endl;
	std::cout << grid.rhs().getSize() << ":" << Amat_.mvmult(SOL).getSize() << std::endl;
	/** Initialize the CG loop **/
	Res = -b - Amat_.mvmult(SOL); // res = b-Ax
	d = Res;
	dnew = Res.dotNC(Res);
	d0 = dnew;
	std::cout << Res.getSize(0) << std::endl;
	std::cout << SOL.getSize(0) << std::endl;
	while(iterno < itermax && resid > eps*p0){
		q = Amat_.mvmult(d);
		alphak = dnew/(d.dotNC(q));
		SOL = SOL + d*alphak;
		if (iterno%checkfrequency == 0){
			Res = -b - Amat_.mvmult(SOL);
		}
		else{
			Res = Res - q*alphak;
		}
		dold = dnew;
		dnew = Res.dotNC(Res);
		betak = dnew/dold;
		d = Res + d*betak;
		resid = sqrt(dnew/nfluid);
		++iterno;
		std::cout<<"Iteration no = "<<iterno<< "\tResidual = "<< resid<<"\n";	
	}
	grid.p() = SOL; //grid.p().reshape(SOL,NxG,NyG);

  if (iterno<itermax) {
      std::cout<<"Solution converged in "<<iterno<< " iterations!\n";
      return true; }
  else
      std::cout<<"\nSolution not converged!\n";
  return false;
}

void PCGSolver::mat_assemble( StaggeredGrid& grid ){
	unsigned int NxG = imax+2;
	unsigned int NyG = jmax+2;
	real dx = grid.dx();
	real dy = grid.dy();
	real idx2 = 1./(dx*dx);
	real idy2 = 1./(dy*dy);
//	MatrixCOO &Amat = Amat_;
	unsigned int index = 0;
	for (unsigned int i = 0; i<NxG; ++i)
		for (unsigned int j = 0; j<NyG; ++j){
			index = i*NyG + j;
			if((i!=0)&&(i!=NxG-1)&&(j!=0)&&(j!=NyG-1)){
				if(grid.isFluid(i,j)){
					Amat_.mat_set(index,index,-2*idx2-2*idy2);
					Amat_.mat_set(index,index+1,idx2);
					Amat_.mat_set(index,index-1,idx2);
					Amat_.mat_set(index,index+NyG,idy2);
					Amat_.mat_set(index,index-NyG,idy2);
				}
			}
			if((i==0)&&(j!=0)&&(j!=NyG-1)){
				Amat_.mat_set(index,index,1);
				Amat_.mat_set(index,index+1,-1);
			}
			if((i==NxG-1)&&(j!=0)&&(j!=NyG-1)){
				Amat_.mat_set(index,index,1);
				Amat_.mat_set(index,index-1,-1);
			}
			if((j==0)&&(i!=0)&&(i!=NxG-1)){
				Amat_.mat_set(index,index,1);
				Amat_.mat_set(index,index+NyG,-1);
			}
			if((j==NyG-1)&&(i!=0)&&(i!=NxG-1)){
				Amat_.mat_set(index,index,1);
				Amat_.mat_set(index,index-NyG,-1);
			}
			if( ((i==0)&&(j==0))||((i==0)&&(j==NyG-1))||((i==NxG-1)&&(j==0))||((i==NxG-1)&&(j==NyG-1)) ){
				Amat_.mat_set(index,index,1);
			}
		}
}
