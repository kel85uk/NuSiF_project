#ifndef STAGGERED_GRID_HH
#define STAGGERED_GRID_HH

#include "Array.hh"
#include "FileReader.hh"
#include "GrayScaleImage.hh"


//*******************************************************************************************************************
/*! Class for storing all Arrays required for the simulation
*
*
* Feel free to add member functions or variables to this class, but please don't
* delete or rename any existing functions, since future skeletons rely on these functions.
*
*/
//*******************************************************************************************************************
typedef enum {NORTH, EAST, SOUTH, WEST, CENTER} Direction;

class StaggeredGrid
{
public:
   // Constructors to manually create staggered grid
   StaggeredGrid ( int xsize, int ysize, real delx, real dely );
   
   // Constructors to create a staggered grid using png file
   StaggeredGrid (GrayScaleImage image, real delx, real dely);

   // Constructor to create a staggered grid from a parsed configuration file
   StaggeredGrid ( const FileReader & config );

   // Getters / Setters for member variables
   Array<real> & p()    { return p_;    }
   Array<real> & u()    { return u_;    }
   Array<real> & v()    { return v_;    }
   Array<real> & f()    { return f_;    }
   Array<real> & g()    { return g_;    }
   Array<real> & rhs()  { return rhs_;  }

   const Array<real> & p()   const { return p_;   }
   const Array<real> & u()   const { return u_;   }
   const Array<real> & v()   const { return v_;   }
   const Array<real> & f()   const { return f_;   }
   const Array<real> & g()   const { return g_;   }
   const Array<real> & rhs() const { return rhs_; }

   real dx() const { return dx_; }
   real dy() const { return dy_; }
   
   const int xSize() const {return p_.getSize(0) - 2;}
   const int ySize() const {return p_.getSize(1) - 2;}
   
   inline bool isFluid(const int x, const int y) {return isfluid_(x,y);}

   inline int getNumFluid();
   inline real& u( int x, int y, Direction dir);
   inline real& v( int x, int y, Direction dir);
   inline real& f( int x, int y, Direction dir);
   inline real& g( int x, int y, Direction dir);
   inline real& p( int x, int y, Direction dir);

protected:
   Array<real> p_;   //< pressure field
   Array<real> u_;   //< velocity field in x direction
   Array<real> v_;   //< velocity field in y direction
   Array<real> f_;   //< 
   Array<real> g_;   //< 
   Array<real> rhs_; //< right hand side of the pressure equation
   Array<unsigned char> isfluid_;

   real xlength,ylength;
   real dx_;   //< distance between two grid points in x direction
   real dy_;   //< distance between two grid points in y direction
   void patch(Array<real> & A, real x1, real y1, real x2, real y2, real value);
   void setCellToObstacle(int x, int y);
   void createRectangle(real x1, real y1, real x2, real y2);
   void createCircle(real x, real y, real r);
};

//===================================================================================================================
//
//  Inline Access Operators and Sizes
//
//===================================================================================================================


// Function to compute the number of fluid cells 
inline int StaggeredGrid::getNumFluid()
{
   int numFluid_ = 0;
   for (int i=1; i<isfluid_.getSize(0)-1; i++) 
       for (int j=1; j<isfluid_.getSize(1)-1; j++)
           numFluid_+=isfluid_(i,j);
   return numFluid_;
}

inline real& StaggeredGrid::u( int x, int y, Direction dir)
{
   switch(dir)
   {
      case NORTH: if (!isfluid_(x,y+1))
                     u_(x,y+1) = -u_(x,y);
                  return u_(x,y+1);
      case SOUTH: if (!isfluid_(x,y-1))
                     u_(x,y-1) = -u_(x,y);
                  return u_(x,y-1);
      case EAST:  if (isfluid_(x+1,y))
                     return u(x+1,y,CENTER);
                  else
                     u_(x+1,y) = 0.0;
                  return u_(x+1,y);
      case WEST:  if (!isfluid_(x-1,y))
                     u_(x-1,y) = 0.0;
                  return u_(x-1,y);
      default:    if (isfluid_(x+1,y))
                     return u_(x,y);
                  else 
                     u_(x,y) = 0.0;
                  return u_(x,y);
   }
}
inline real& StaggeredGrid::v(int x, int y, Direction dir)
{
   switch(dir)
   {
      case NORTH: if (isfluid_(x,y+1))
                     return v(x,y+1,CENTER);
                  else
                     v_(x,y+1)=0.0;
                     return v_(x,y+1);
      case SOUTH: if (!isfluid_(x,y-1))
                     v_(x,y-1)=0.0;
                  return v_(x,y-1);
      case EAST:  if (!isfluid_(x+1,y))
                     v_(x+1,y)=-v_(x,y);
                  return v_(x+1,y);
      case WEST:  if (isfluid_(x-1,y))
                     v_(x-1,y)=-v_(x,y);
                  return v_(x-1,y);
      default:    if (!isfluid_(x,y+1))
                     v_(x,y) = 0.0;
                  return v_(x,y);
   }
}

inline real& StaggeredGrid::f( int x, int y, Direction dir)
{
   switch(dir)
   {
      case CENTER: if (isfluid_(x+1,y))
                     return f_(x,y);
                  else
                     return u_(x,y);
      case WEST:  if (isfluid_(x-1,y))
                     return f_(x-1,y);
                  else
                     return u_(x-1,y);
      default:    return f_(x,y);
   }
}
inline real& StaggeredGrid::g(int x, int y, Direction dir)
{
   switch(dir)
   {
      case CENTER: if (isfluid_(x,y+1))
                     return g_(x,y);
                  else
                     return v_(x,y);
      case SOUTH: if (isfluid_(x,y-1))
                     return g_(x,y-1);
                  else
                     return v_(x,y-1);
      default:    return g_(x,y);
   }
}

inline real& StaggeredGrid::p(int x, int y, Direction dir)
{
   switch(dir)
   {
      case NORTH: if (isfluid_(x,y+1))
                     return p_(x,y+1);
                  else
                     return p_(x,y);
      case SOUTH: if (isfluid_(x,y-1))
                     return p_(x,y-1);
                  else
                     return p_(x,y);
      case EAST:  if (isfluid_(x+1,y))
                     return p_(x+1,y);
                  else
                     return p_(x,y);
      case WEST:  if (isfluid_(x-1,y))
                     return p_(x-1,y);
                  else
                     return p_(x,y);
      default:    return p_(x,y);
   }
}

#endif //STAGGERED_GRID_HH
