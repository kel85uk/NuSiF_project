#include "StaggeredGrid3D.hh"

//===================================================================================================================
//
//  Constructors
//
//===================================================================================================================

   // Constructors to manually create staggered grid
StaggeredGrid3D::StaggeredGrid3D (int xsize, int ysize, int zsize, real delx, real dely, real delz):
					p_(xsize+2,ysize+2,zsize+2), u_(xsize+1,ysize+2,zsize+2), v_(xsize+2,ysize+1,zsize+2), 
					w_(xsize+2,ysize+1,zsize+2), f_(xsize+1,ysize+1,zsize+1), g_(xsize+1,ysize+1,zsize+1),
					h_(xsize+1,ysize+1,zsize+1), rhs_(xsize+2,ysize+2,zsize+2),
					isfluid_(xsize+2,ysize+2,zsize+2), dx_(delx), dy_(dely), dz_(delz)
{
   std::cout<<"\nMaking Staggered Grid\n";
   p_.fill(1.0);
   u_.fill(0.0);
   v_.fill(0.0);
   w_.fill(0.0);
   isfluid_.fill(255);
}

   // Constructor to create a staggered grid from a parsed config file
StaggeredGrid3D::StaggeredGrid3D (const FileReader & config): 
p_(config.getIntParameter("imax")+2, config.getIntParameter("jmax")+2, config.getIntParameter("kmax")+2),
u_(config.getIntParameter("imax")+1, config.getIntParameter("jmax")+2, config.getIntParameter("kmax")+2),
v_(config.getIntParameter("imax")+2, config.getIntParameter("jmax")+1, config.getIntParameter("kmax")+2),
w_(config.getIntParameter("imax")+2, config.getIntParameter("jmax")+2, config.getIntParameter("kmax")+1),
f_(config.getIntParameter("imax")+1, config.getIntParameter("jmax")+1, config.getIntParameter("kmax")+1),
g_(config.getIntParameter("imax")+1, config.getIntParameter("jmax")+1, config.getIntParameter("kmax")+1),
h_(config.getIntParameter("imax")+1, config.getIntParameter("jmax")+1, config.getIntParameter("kmax")+1),
rhs_(config.getIntParameter("imax")+1, config.getIntParameter("jmax")+1, config.getIntParameter("kmax")+1),
isfluid_(config.getIntParameter("imax")+2, config.getIntParameter("jmax")+2, config.getIntParameter("kmax")+2),
xlength(config.getRealParameter("xlength")), ylength(config.getRealParameter("ylength")),
zlength(config.getRealParameter("zlength")), dx_(xlength/config.getIntParameter("imax")),
dy_(ylength/config.getIntParameter("jmax")), dz_(zlength/config.getIntParameter("kmax"))
{
   std::cout<<"\nMaking Staggered Grid\n";
   isfluid_.fill(255);
   for (int i=2; i<isfluid_.getSize(0)-1; i++) 
       for (int j=2; j<isfluid_.getSize(1)-1; j++)
           for (int k=2; k<isfluid_.getSize(2)-1; k++)
               if (!isfluid_(i,j,k))
                   CHECK_MSG(!((isfluid_(i,j+1,k) && isfluid_(i,j-1,k)) || (isfluid_(i+1,j,k) && isfluid_(i-1,j,k))
                                 || (isfluid_(i,j,k+1) && isfluid_(i,j,k-1))), "Walls must be minimum 2 cell thick");

   std::cout<<"\nInitializing Grid\n";

   if (!config.find("P_INIT")) {
       WARN("\nInitial value of P not specified\n");
       p_.fill(1.0); }
   else 
       p_.fill(config.getRealParameter("P_INIT"));
       
   if (!config.find("U_INIT")) {
       WARN("\nInitial value of U not specified\n");
       u_.fill(0.0); }
   else 
       u_.fill(config.getRealParameter("U_INIT"));
       
   if (!config.find("V_INIT")) {
       WARN("\nInitial value of V not specified\n");
       v_.fill(0.0); }
   else 
       v_.fill(config.getRealParameter("V_INIT"));

   if (!config.find("W_INIT")) {
       WARN("\nInitial value of W not specified\n");
       w_.fill(0.0); }
   else 
       w_.fill(config.getRealParameter("W_INIT"));
}

void StaggeredGrid3D::setCellToObstacle(int x, int y, int z)
{
   isfluid_(x,y,z) = 0;
}
