#include "StaggeredGrid.hh"

//===================================================================================================================
//
//  Constructors
//
//===================================================================================================================

   // Constructors to manually create staggered grid
StaggeredGrid::StaggeredGrid (int xsize, int ysize, real delx, real dely):
					p_(xsize+2,ysize+2), u_(xsize+1,ysize+2), v_(xsize+2,ysize+1),
					f_(xsize+1,ysize+1), g_(xsize+1,ysize+1), rhs_(xsize+2,ysize+2),
					isfluid_(xsize+2,ysize+2), dx_(delx), dy_(dely)
{
   std::cout<<"\nMaking Staggered Grid\n";
   isfluid_.fill(255);
   std::cout<<"\nInitializing Grid\n";
   p_.fill(1.0);
   u_.fill(0.0);
   v_.fill(0.0);
   rhs_.fill(0.0);
}

   // Constructors to create staggered grid from png file
StaggeredGrid::StaggeredGrid (GrayScaleImage image, real delx, real dely):
					p_(image.width()+2,image.height()+2), u_(image.width()+1,image.height()+2),
					v_(image.width()+2,image.height()+1), f_(image.width()+1,image.height()+1), 
					g_(image.width()+1,image.height()+1), rhs_(image.width()+2,image.height()+2),
					isfluid_(image.width()+2,image.height()+2), dx_(delx), dy_(dely)
{
   std::cout<<"\nMaking Staggered Grid\n";
   std::cout<<"Reading Image file\n";
   for (int i=1; i<isfluid_.getSize(0); i++) 
       for (int j=1; j<isfluid_.getSize(1); j++)
           isfluid_(i,j)= image.getElement(i,j);
   for (int i=2; i<isfluid_.getSize(0)-1; i++) 
       for (int j=2; j<isfluid_.getSize(1)-1; j++)
           if (!isfluid_(i,j))
               CHECK_MSG(!((isfluid_(i,j+1) && isfluid_(i,j-1)) || (isfluid_(i+1,j) && isfluid_(i-1,j))),
                         "Walls must be minimum 2 cell thick") ;

   // Counting cells
   numFluid_ = 0;
   for (int i=1; i<isfluid_.getSize(0)-1; i++) 
       for (int j=1; j<isfluid_.getSize(1)-1; j++)
           numFluid_+=isfluid_(i,j);

   std::cout<<"\nInitializing Grid\n";
   p_.fill(1.0);
   u_.fill(0.0);
   v_.fill(0.0);
}

   // Constructor to create a staggered grid from a parsed configuration file
StaggeredGrid::StaggeredGrid (const FileReader & config): 
                              p_(config.getIntParameter("imax")+2, config.getIntParameter("jmax")+2),
                              u_(config.getIntParameter("imax")+1, config.getIntParameter("jmax")+2),
                              v_(config.getIntParameter("imax")+2, config.getIntParameter("jmax")+1),
                              f_(config.getIntParameter("imax")+1, config.getIntParameter("jmax")+1),
                              g_(config.getIntParameter("imax")+1, config.getIntParameter("jmax")+1),
                              rhs_(config.getIntParameter("imax")+2 ,config.getIntParameter("jmax")+2),
                              isfluid_(config.getIntParameter("imax")+2 ,config.getIntParameter("jmax")+2),
                              xlength(config.getRealParameter("xlength")), ylength(config.getRealParameter("ylength")),
                              dx_(xlength/config.getIntParameter("imax")), dy_(ylength/config.getIntParameter("jmax"))

{
   std::cout<<"\nMaking Staggered Grid\n";
   isfluid_.fill(255);
   if (config.find("RectangleX1"))
   {
       createRectangle(config.getRealParameter("RectangleX1"), config.getRealParameter("RectangleY1"),
                       config.getRealParameter("RectangleX2"), config.getRealParameter("RectangleY2"));
   }
   if (config.find("CircleX"))
   {
       createCircle(config.getRealParameter("CircleX"), config.getRealParameter("CircleY"),
                    config.getRealParameter("CircleR")); 
   }
   for (int i=2; i<isfluid_.getSize(0)-1; i++) 
       for (int j=2; j<isfluid_.getSize(1)-1; j++)
           if (!isfluid_(i,j))
               CHECK_MSG(!((isfluid_(i,j+1) && isfluid_(i,j-1)) || (isfluid_(i+1,j) && isfluid_(i-1,j))),
                         "Walls must be minimum 2 cell thick") ;
 
   // Counting cells
   numFluid_ = 0;
   for (int i=1; i<isfluid_.getSize(0)-1; i++) 
       for (int j=1; j<isfluid_.getSize(1)-1; j++)
           numFluid_+=isfluid_(i,j);

   std::cout<<"\nWriting Grid to file\n";
   
   GrayScaleImage image(isfluid_);
   image.save(config.getStringParameter("name"));

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
   std::cout<<"\nPatching Grid\n";
   if (config.find("UpatchX1"))
   {
       patch(u_, config.getRealParameter("UpatchX1"), config.getRealParameter("UpatchY1"),
             config.getRealParameter("UpatchX2"), config.getRealParameter("UpatchY2"),
             config.getRealParameter("UpatchValue"));
   }
   if (config.find("VpatchX1"))
   {
       patch(v_, config.getRealParameter("VpatchX1"), config.getRealParameter("VpatchY1"),
             config.getRealParameter("VpatchX2"), config.getRealParameter("VpatchY2"),
             config.getRealParameter("VpatchValue"));
   }
   if (config.find("PpatchX1"))
   {
       patch(p_, config.getRealParameter("PpatchX1"), config.getRealParameter("PpatchY1"),
             config.getRealParameter("PpatchX2"), config.getRealParameter("PpatchY2"),
             config.getRealParameter("PpatchValue"));
   }
}

void StaggeredGrid::patch(Array<real> & A, real x1, real y1, real x2, real y2, real value)
{
   CHECK_MSG(x1 >= 0 && x2 > x1 && y1>=0 && y2>y1,
             "Patch geometry must be specified by setting the bottom-left and top-right point");
   CHECK_MSG(x2 <= xlength && y2 <= ylength, "Patch must lie inside the domain");
   for (int i=1; i<p_.getSize(0)-1; i++) 
       for (int j=1; j<p_.getSize(1)-1; j++)
           if(dx_*(i-0.5) >= x1 && dx_*(i-0.5) <= x2 && dy_*(j-0.5) >= y1 && dy_*(j-0.5) <= y2)
              A(i,j) = value;
}   

void StaggeredGrid::setCellToObstacle(int x, int y)
{
   isfluid_(x,y) = 0;
}

void StaggeredGrid::createRectangle(real x1, real y1, real x2, real y2)
{
    CHECK_MSG(x1 >= 0 && x2 > x1 && y1>=0 && y2>y1,
              "Rectangle geometry must be specified by setting the bottom-left and top-right point");
    CHECK_MSG(x2 <= xlength && y2 <= ylength, "Rectangle must lie inside the domain");
   for (int i=1; i<isfluid_.getSize(0); i++) 
       for (int j=1; j<isfluid_.getSize(1); j++)
           if(dx_*(i-0.5) >= x1 && dx_*(i-0.5) <= x2 && dy_*(j-0.5) >= y1 && dy_*(j-0.5) <= y2)
              setCellToObstacle(i,j);
}   

void StaggeredGrid::createCircle(real x, real y, real r)
{
   for (int i=1; i<isfluid_.getSize(0)-1; i++) 
       for (int j=1; j<isfluid_.getSize(1)-1; j++)
           if( pow((x-(i-0.5)*dx_),2) + pow(y-j*dy_,2) <= r*r)
              setCellToObstacle(i,j);
}
