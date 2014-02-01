#include "MultigridSolver.hh"

MultigridSolver::MultigridSolver ( unsigned int imax_, unsigned int jmax_, real omg_, real eps_, unsigned int itermax_):
                     imax(imax_), jmax(jmax_), itermax(itermax_), checkfrequency(1), omg(omg_), eps(eps_)
{
   std::cout<<"\nCreating Multigrid solver\n";
}

MultigridSolver::MultigridSolver ( const FileReader & configuration ):
                     imax(configuration.getIntParameter("imax")), jmax(configuration.getIntParameter("jmax")),
                     itermax(configuration.getIntParameter("itermax")),
                     levels(configuration.getIntParameter("MultigridLevels")), itercoarse(configuration.getIntParameter("itercoarse")),
                     iterpre(configuration.getIntParameter("iterpre")),iterpost(configuration.getIntParameter("iterpost")),
                     omg(configuration.getRealParameter("omg")), eps(configuration.getRealParameter("eps"))
{
   std::cout<<"\nCreating Multigrid solver\n";
   CHECK_MSG((omg <= 1.9 && omg >= 0.5), "omg must be between 1.7 and 1.9");
   CHECK_MSG((eps > 0), "eps must be greater than 0");
   if (!configuration.find ("checkfrequency")) {
       WARN("\nCheck Frequency not specified\n");
      checkfrequency = 1; }
   else  
      checkfrequency = configuration.getIntParameter("checkfrequency");
   CHECK_MSG((checkfrequency > 0), "check frequency must be greater than 0");
   CHECK_MSG((itermax > 0), "Max number of iterations (itermax) must be greater than 0");
}

void MultigridSolver::SetBoundary ( Array<real> & p )
{
   for (int i=0; i < p.getSize(0); i++) {
        p(i,0) = p(i,1);
        p(i, p.getSize(1) - 1) = p(i, p.getSize(1) - 2); 
   }
   for (int j=0; j< p.getSize(1); j++) {
        p(0,j) = p(1,j);
        p(p.getSize(0) - 1, j) = p(p.getSize(0) - 2, j);
   }
}

// Calculate residual
inline real MultigridSolver::Residual ( Array<real> & sol, Array<real> & rhs, real dx, real dy )
{
   real resid=0.0, dx_2=1.0/(dx*dx), dy_2=1.0/(dy*dy), res;

   for (int i=1; i<sol.getSize(0)-1; i++){
       for (int j=1; j<sol.getSize(1)-1; j++){
	   	res = (sol(i+1,j) + sol(i-1,j) - 2.0*sol(i,j))*dx_2 + (sol(i,j+1) + sol(i,j-1) - 2.0*sol(i,j))*dy_2 - rhs(i,j) ;
        resid += res*res;
	}
   }

   return ( sqrt( resid/(sol.getSize(0)*sol.getSize(1)) ) );
}


bool MultigridSolver::solve( StaggeredGrid & grid )
{
   	Array<real> & p = grid.p();
   	Array<real> & rhs = grid.rhs();
   	real res , dx = grid.dx(), dy = grid.dy();
   	int nrows = grid.xSize(), ncols = grid.ySize();
   
	Array <real>** Sol = new Array<real>*[levels];
	Array <real>** RHS = new Array<real>*[levels];

	SetBoundary(p);
   
	Sol[0] = new Array<real>(p);
	RHS[0] = new Array<real>(rhs);
		
	for (unsigned int i = 1; i < levels; i++ ){
		nrows = nrows/2 ;
		ncols = ncols/2 ;
		
		Sol[i] = new Array<real>( nrows+2 ,ncols+2 );
		RHS[i] = new Array<real>( nrows+2 ,ncols+2 );
	}
	
	for ( unsigned int i = 1; i < levels; i++ ){
		Restrict( *Sol[i-1], *Sol[i] );
		Restrict( *RHS[i-1], *RHS[i] );
	}
	std::cout<<"\nStarting Multigrid Iteration.\n";  
	// Multigrid iteration
	unsigned int iterno = 0;
   //	for ( unsigned int lvl = levels; lvl != 0; lvl-- ) {
   //		unsigned int lev = lvl-1; 
   		unsigned int lev = 0;
   		
   		for ( iterno = 0; iterno < itermax; iterno++ ){ 			
			VCycle ( lev,Sol,RHS,dx,dy);
	
			if (iterno%checkfrequency == 0 ) {
				res = Residual ( *Sol[lev], *RHS[lev] , pow(2,lev)*dx, pow(2,lev)*dy);
				std::cout<<"Level = " << lev+1 <<"\tIteration no = "<<iterno<< "\tResidual = "<< res<<"\n"; 
				if (res < eps) break;
			}
		}

	//	if ( lev > 0 )
	//		interpolate ( *Sol[lev - 1],*Sol[lev]);	
	//}

   	for (int i=0; i<nrows+2; i++){
   		for (int j=0; j<ncols+2; j++){
   			
   			p(i,j) = (*Sol[0])(i,j) ;
	   		rhs(i,j) = (*RHS[0])(i,j) ;
		}
   	}
   	
   	for (unsigned int i = 0; i < levels; i++ ){
		
		delete Sol[i];
		delete RHS[i];
	}

	delete[] Sol;
	delete[] RHS;

  	if (iterno<itermax) {
    	std::cout<<"Solution converged in "<<iterno<< " iterations!\n";
      	return true;
    }
  	else{
      	std::cout<<"\nSolution not converged!\n";
  		return false;
  	}
}

void MultigridSolver::VCycle(unsigned int lev, Array<real>** sol, Array<real>** rhs, real dx, real dy){
	real local_dx = pow(2,lev)*dx, local_dy = pow(2,lev)*dx; 

	if(lev == levels-1){
		for(unsigned int i=0; i < itercoarse; i++){
			
			SORsweep(*sol[lev], *rhs[lev], local_dx, local_dy);
			real res = Residual ( *sol[lev], *rhs[lev] , local_dx, local_dy);
			//std::cout<<"Internal : Iteration no = "<<i<< "\tResidual = "<< res<<"\n"; 
			if (res < eps) break;
		}
	}
	else {
		for(unsigned int i=0; i < iterpre; i++)
			SORsweep(*sol[lev], *rhs[lev], local_dx, local_dy);

		RestrictResisdual(*sol[lev],*rhs[lev],*rhs[lev+1], local_dx, local_dy);
			
//std::cout<<"RHS summary = " << rhsSum(*rhs[lev+1]) <<std::endl;		
		
		sol[lev+1]->fill(0);
		
		VCycle(lev+1, sol, rhs, dx, dy);
		
		interpolateCorrect(*sol[lev], *sol[lev+1]);
		
		for(unsigned int i=0; i < iterpost; i++)
			SORsweep(*sol[lev], *rhs[lev], local_dx, local_dy);
	}
}

real MultigridSolver::rhsSum( Array<real> & rhs ){
	
   	real sum = 0.0;
   	for (int i=1; i<rhs.getSize(0)-1; i++)
    	for (int j=1; j<rhs.getSize(1)-1; j++)
        	sum += rhs(i,j);
               
	return sum;
}

void MultigridSolver::SORsweep( Array<real> & sol, Array<real> & rhs, real dx, real dy ){
	real dx_2 = 1.0/(dx*dx), dy_2 = 1.0/(dy*dy);
	real facD = omg / (2.0/(dx*dx) + 2.0/(dy*dy)) ;
	SetBoundary(sol);
	
	for (int i=1; i<sol.getSize(0)-1; i++)
        	for (int j=1; j<sol.getSize(1)-1; j++)
                	sol(i,j) = (1.0-omg)* sol(i,j) 
					+ facD * ( (sol(i+1,j) + sol(i-1,j) )*dx_2 + ( sol(i,j+1) + sol(i,j-1) )*dy_2 - rhs(i,j) );
    SetBoundary(sol);

}

void MultigridSolver::Restrict (Array<real>& fine, Array<real>& coarse){

	for (int i=1; i<coarse.getSize(0)-1 ;i++){
		int fi = 2*i;
		for (int j=1; j<coarse.getSize(1)-1 ;j++){
			int fj = 2*j;

			coarse(i,j) = 0.25 * ( fine(fi,fj) + fine(fi-1,fj) + fine(fi-1,fj-1) + fine(fi,fj-1) );
		}
	}
}

void MultigridSolver::RestrictResisdual(Array<real>& sol, Array<real>& rhs, Array<real>& rhs_coarse, real dx, real dy){
	
	real dx_2=1.0/(dx*dx), dy_2=1.0/(dy*dy);

	for (int i=1; i<rhs_coarse.getSize(0)-1 ;i++){
		int fi = 2*i;
		for (int j=1; j<rhs_coarse.getSize(1)-1 ;j++){
			int fj = 2*j;
			
			real res_ij = (sol(fi+1,fj) + sol(fi-1,fj) - 2.0*sol(fi,fj))*dx_2 + (sol(fi,fj+1) + sol(fi,fj-1) - 2.0*sol(fi,fj))*dy_2 - rhs(fi,fj) ;
			real res_i1j = (sol(fi,fj) + sol(fi-2,fj) - 2.0*sol(fi-1,fj))*dx_2 + (sol(fi-1,fj+1) + sol(fi-1,fj-1) - 2.0*sol(fi-1,fj))*dy_2 - rhs(fi-1,fj) ;
			real res_ij1 = (sol(fi+1,fj-1) + sol(fi-1,fj-1) - 2.0*sol(fi,fj-1))*dx_2 + (sol(fi,fj) + sol(fi,fj-2) - 2.0*sol(fi,fj-1))*dy_2 - rhs(fi,fj-1);
			real res_i1j1 = (sol(fi,fj-1) + sol(fi-2,fj-1) - 2.0*sol(fi-1,fj-1))*dx_2 + (sol(fi-1,fj) + sol(fi-1,fj-2) - 2.0*sol(fi-1,fj-1))*dy_2 - rhs(fi-1,fj-1);

			rhs_coarse(i,j) = 0.25 * ( res_ij + res_i1j + res_ij1 + res_i1j1 );
		}
	}
}

void MultigridSolver::interpolateCorrect(Array<real>& sol_fine, Array<real>& sol_coarse){

	double v = 0;

	for (int i=1; i<sol_coarse.getSize(0)-1 ;i++){
		int fi = 2*i;
		for (int j=1; j<sol_coarse.getSize(1)-1 ;j++){
			int fj = 2*j;
			v = sol_coarse(i,j);
				
			sol_fine( fi,fj ) += v;
			sol_fine( fi-1,fj ) += 0.5*( sol_coarse(i+1,j) + sol_coarse(i-1,j) );
			sol_fine( fi,fj-1 ) += 0.5*( sol_coarse(i,j+1) + sol_coarse(i,j-1) );
			sol_fine( fi-1,fj-1 ) += 0.25*( sol_coarse(i+1,j) + sol_coarse(i-1,j) + sol_coarse(i,j+1) + sol_coarse(i,j-1) );
		}
	}
}

void MultigridSolver::interpolate(Array<real>& sol_fine, Array<real>& sol_coarse){
	
	double v = 0;
	
	for (int i=1; i<sol_coarse.getSize(0)-1 ;i++){
		int fi = 2*i;
		for (int j=1; j<sol_coarse.getSize(1)-1 ;j++){
			int fj = 2*j;

			v = sol_coarse(i,j);
			
			sol_fine( fi,fj ) = v;
			sol_fine( fi-1,fj ) = 0.5*( sol_coarse(i+1,j) + sol_coarse(i-1,j) );
			sol_fine( fi,fj-1 ) = 0.5*( sol_coarse(i,j+1) + sol_coarse(i,j-1) );
			sol_fine( fi-1,fj-1 ) = 0.25*( sol_coarse(i+1,j) + sol_coarse(i-1,j) + sol_coarse(i,j+1) + sol_coarse(i,j-1) );
		}
	}
}
