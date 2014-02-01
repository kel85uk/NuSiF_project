#ifndef STAGGERED_GRID_3D_HH
#define STAGGERED_GRID_3D_HH

#include "Array.hh"
#include "FileReader.hh"


//*******************************************************************************************************************
/*! Class for storing all arrays required for the simulation
*
* For now it only contains an array for the pressure and another array for the
* right hand side of the pressure equation.
* In following assignments this will be extended and also contain
* arrays for x/y velocity and some arrays for intermediate values.
*
* Feel free to add member functions or variables to this class, but please don't
* delete or rename any existing functions, since future skeletons rely on these functions.
*
*/
//*******************************************************************************************************************
typedef enum {NORTH, EAST, SOUTH, WEST, UP, DOWN, CENTER} Direction;

class StaggeredGrid3D
{
public:
   // Constructors to manually create staggered grid
   StaggeredGrid3D ( int xsize, int ysize, int zsize, real delx, real dely, real delz);

   // Constructor to create a staggered grid from a parsed configuration file
   StaggeredGrid3D ( const FileReader & config );

   // Getters / Setters for member variables
   Array<real> & p()    { return p_;    }
   Array<real> & u()    { return u_;    }
   Array<real> & v()    { return v_;    }
   Array<real> & w()    { return w_;    }
   Array<real> & f()    { return f_;    }
   Array<real> & g()    { return g_;    }
   Array<real> & h()    { return h_;    }
   Array<real> & rhs()  { return rhs_;  }

   const Array<real> & p()   const { return p_;   }
   const Array<real> & u()   const { return u_;   }
   const Array<real> & v()   const { return v_;   }
   const Array<real> & w()   const { return w_;   }
   const Array<real> & f()   const { return f_;   }
   const Array<real> & g()   const { return g_;   }
   const Array<real> & h()   const { return h_;   }
   const Array<real> & rhs() const { return rhs_; }

   real dx() const { return dx_; }
   real dy() const { return dy_; }
   real dz() const { return dz_; }
   
   const int xSize() const {return p_.getSize(0) - 2;}
   const int ySize() const {return p_.getSize(1) - 2;}
   const int zSize() const {return p_.getSize(2) - 2;}

   inline bool isFluid(const int x, const int y, const int z) {return isfluid_(x,y,z);}

   inline int getNumFluid();
   inline real& u( int x, int y, int z, Direction dir);
   inline real& v( int x, int y, int z, Direction dir);
   inline real& w( int x, int y, int z, Direction dir);
   inline real& f( int x, int y, int z, Direction dir);
   inline real& g( int x, int y, int z, Direction dir);
   inline real& h( int x, int y, int z, Direction dir);
   inline real& p( int x, int y, int z, Direction dir);


protected:
   Array<real> p_;   //< pressure field
   Array<real> u_;   //< velocity field in x direction
   Array<real> v_;   //< velocity field in y direction
   Array<real> w_;   //< velocity field in z direction
   Array<real> f_;   //< 
   Array<real> g_;   //< 
   Array<real> h_;   //< 
   Array<real> rhs_; //< right hand side of the pressure equation
   Array<unsigned char> isfluid_;

   int numFluid_;
   real xlength,ylength,zlength;

   real dx_;   //< distance between two grid points in x direction
   real dy_;   //< distance between two grid points in y direction
   real dz_;   //< distance between two grid points in z direction
   
   void setCellToObstacle(int x, int y, int z);
};


//===================================================================================================================
//
//  Inline Access Operators and Sizes
//
//===================================================================================================================

// Function to compute the number of fluid cells 
inline int StaggeredGrid3D::getNumFluid()
{
   return numFluid_;
}

inline real& StaggeredGrid3D::u( int x, int y, int z, Direction dir)
{
   switch(dir)
   {
      case NORTH: if (!(isfluid_(x,y+1,z) || isfluid_(x+1,y+1,z)))
                      u_(x,y+1,z) = -u_(x,y,z);
                  else if (isfluid_(x,y+1,z) && isfluid_(x+1,y+1,z))
                      return u_(x,y+1,z);
                  else
                      u_(x,y+1,z) = 0;
                  return u_(x,y+1);
      case SOUTH:  if (!(isfluid_(x,y-1,z) || isfluid_(x+1,y-1,z)))
                      u_(x,y-1,z) = -u_(x,y,z);
                  else if (isfluid_(x,y-1,z) && isfluid_(x+1,y-1,z))
                      return u_(x,y-1,z);
                  else
                      u_(x,y-1,z) = 0;
                  return u_(x,y+1,z);
      case EAST:  if (isfluid_(x+1,y,z))
                     return u(x+1,y,z,CENTER);
                  else
                     u_(x+1,y,z) = 0.0;
                     return u_(x+1,y,z);
      case WEST:  if (isfluid_(x-1,y,z))
                     return u_(x-1,y,z);
                  else
                     u_(x-1,y,z) = 0.0;
                     return u_(x-1,y,z);
      case UP: if (!(isfluid_(x,y,z+1) || isfluid_(x,y+1,z+1)))
                      u_(x,y,z+1) = -u_(x,y,z);
                  else if (isfluid_(x,y,z+1) && isfluid_(x,y+1,z+1))
                      return u_(x,y,z+1);
                  else
                      u_(x,y,z+1) = 0;
                  return u_(x,y,z+1);
      case DOWN: if (!(isfluid_(x,y,z-1) || isfluid_(x+1,y,z-1)))
                      u_(x,y,z-1) = -u_(x,y,z);
                  else if (isfluid_(x,y,z-1) && isfluid_(x+1,y,z-1))
                      return u_(x,y,z-1);
                  else
                      u_(x,y,z-1) = 0;
                  return u_(x,y,z-1);
      default:    if (isfluid_(x+1,y,z))
                     return u_(x,y,z);
                  else 
                     u_(x,y,z) = 0.0;
                     return u_(x,y,z);
   }
}

inline real& StaggeredGrid3D::v(int x, int y, int z, Direction dir)
{
   switch(dir)
   {
      case NORTH: if (isfluid_(x,y+1,z))
                     return v(x,y+1,z,CENTER);
                  else
                     v_(x,y+1,z)=0;
                     return v_(x,y+1,z);
      case SOUTH: if (isfluid_(x,y-1,z))
                     return v_(x,y-1,z);
                  else
                     v_(x,y-1,z) = 0.0;
                     return v_(x,y-1,z);
      case EAST: if (!(isfluid_(x+1,y,z) || isfluid_(x+1,y+1,z)))
                      v_(x+1,y,z) = -v_(x,y,z);
                  else if (isfluid_(x+1,y,z) && isfluid_(x+1,y+1,z))
                      return v_(x+1,y,z);
                  else
                      v_(x+1,y,z) = 0;
                  return v_(x+1,y,z);
      case WEST: if (!(isfluid_(x-1,y,z) || isfluid_(x-1,y+1,z)))
                      v_(x-1,y,z) = -v_(x,y,z);
                  else if (isfluid_(x-1,y,z) && isfluid_(x-1,y+1,z))
                      return v_(x-1,y,z);
                  else
                      v_(x-1,y,z) = 0;
                  return v_(x-1,y,z);
      case UP: if (!(isfluid_(x,y,z+1) || isfluid_(x,y+1,z+1)))
                      v_(x,y,z+1) = -v_(x,y,z);
                  else if (isfluid_(x,y,z+1) && isfluid_(x,y+1,z+1))
                      return v_(x,y,z+1);
                  else
                      v_(x,y,z+1) = 0;
                  return v_(x,y,z+1);
      case DOWN: if (!(isfluid_(x,y,z-1) || isfluid_(x,y+1,z-1)))
                      v_(x,y,z-1) = -v_(x,y,z);
                  else if (isfluid_(x,y,z-1) && isfluid_(x,y+1,z-1))
                      return v_(x,y,z-1);
                  else
                      v_(x,y,z-1) = 0;
                  return v_(x,y,z-1);
      default:    if (isfluid_(x,y+1,z))
                     return v_(x,y,z);
                  else 
                     v_(x,y,z) = 0.0;
                     return v_(x,y,z);
   }
}

inline real& StaggeredGrid3D::w(int x, int y, int z, Direction dir)
{
   switch(dir)
   {
      case NORTH: if (!(isfluid_(x,y+1,z) || isfluid_(x,y+1,z+1)))
                      w_(x,y+1,z) = -w_(x,y,z);
                  else if (isfluid_(x,y+1,z) && isfluid_(x,y+1,z+1))
                      return w_(x,y+1,z);
                  else
                      w_(x,y+1,z) = 0;
                  return w_(x,y+1);
      case SOUTH:  if (!(isfluid_(x,y-1,z) || isfluid_(x,y-1,z+1)))
                      w_(x,y-1,z) = -w_(x,y,z);
                  else if (isfluid_(x,y-1,z) && isfluid_(x,y-1,z+1))
                      return w_(x,y-1,z);
                  else
                      w_(x,y-1,z) = 0;
                  return w_(x,y+1,z);
      case EAST: if (!(isfluid_(x+1,y,z) || isfluid_(x+1,y,z+1)))
                      w_(x+1,y,z) = -w_(x,y,z);
                  else if (isfluid_(x+1,y,z) && isfluid_(x+1,y,z+1))
                      return w_(x+1,y,z);
                  else
                      w_(x+1,y,z) = 0;
                  return w_(x+1,y,z);
      case WEST: if (!(isfluid_(x-1,y,z) || isfluid_(x-1,y,z+1)))
                      w_(x-1,y,z) = -w_(x,y,z);
                  else if (isfluid_(x-1,y,z) && isfluid_(x-1,y,z+1))
                      return w_(x-1,y,z);
                  else
                      w_(x-1,y,z) = 0;
                  return w_(x-1,y,z);
      case UP:  if (isfluid_(x,y,z+1))
                     return w(x,y,z+1,CENTER);
                  else
                     w_(x,y,z+1) = 0.0;
                     return w_(x,y,z+1);
      case DOWN:  if (isfluid_(x,y,z-1))
                     return w_(x,y,z-1);
                  else
                     w_(x,y,z-1) = 0.0;
                     return w_(x,y,z-1);
      default:    if (isfluid_(x,y,z+1))
                     return w_(x,y,z);
                  else
                     w_(x,y,z) = 0.0;
                     return w_(x,y,z);
   }
}

inline real& StaggeredGrid3D::f(int x, int y, int z, Direction dir)
{
   switch(dir)
   {
      case CENTER: if (isfluid_(x+1,y,z))
                     return f_(x,y,z);
                  else
                     return u_(x,y,z);
      case WEST:  if (isfluid_(x-1,y,z))
                     return f_(x-1,y,z);
                  else
                     return u_(x-1,y,z);
      default:    return f_(x,y,z);
   }
}

inline real& StaggeredGrid3D::g(int x, int y, int z, Direction dir)
{
   switch(dir)
   {
      case CENTER: if (isfluid_(x,y+1,z))
                     return g_(x,y,z);
                  else
                     return v_(x,y,z);
      case SOUTH: if (isfluid_(x,y-1,z))
                     return g_(x,y-1,z);
                  else
                     return v_(x,y-1,z);
      default:    return g_(x,y,z);
   }
}

inline real& StaggeredGrid3D::h(int x, int y, int z, Direction dir)
{
   switch(dir)
   {
      case CENTER: if (isfluid_(x,y,z+1))
                     return h_(x,y,z);
                  else
                     return w_(x,y,z);
      case DOWN: if (isfluid_(x,y,z-1))
                     return h_(x,y,z-1);
                  else
                     return w_(x,y,z-1);
      default:    return h_(x,y,z);
   }
}

inline real& StaggeredGrid3D::p(int x, int y, int z, Direction dir)
{
   switch(dir)
   {
      case NORTH: if (isfluid_(x,y+1,z))
                     return p_(x,y+1,z);
                  else
                     return p_(x,y,z);
      case SOUTH: if (isfluid_(x,y-1,z))
                     return p_(x,y-1,z);
                  else
                     return p_(x,y,z);
      case EAST:  if (isfluid_(x+1,y,z))
                     return p_(x+1,y,z);
                  else
                     return p_(x,y,z);
      case WEST:  if (isfluid_(x-1,y,z))
                     return p_(x-1,y,z);
                  else
                     return p_(x,y,z);
      case UP:  if (isfluid_(x,y,z+1))
                     return p_(x,y,z+1);
                  else
                     return p_(x,y,z);
      case DOWN:  if (isfluid_(x,y,z-1))
                     return p_(x,y,z-1);
                  else
                     return p_(x,y,z);
      default:    return p_(x,y,z);
   }
}

#endif //STAGGERED_GRID_3D_HH
