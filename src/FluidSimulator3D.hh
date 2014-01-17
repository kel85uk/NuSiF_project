#ifndef __FLUID_SIMULATOR_3D_H__
#define __FLUID_SIMULATOR_3D_H__


#include "FileReader.hh"
#include "StaggeredGrid3D.hh"
#include "SORSolver3D.hh"
class FluidSimulator3D
{
  public:
      FluidSimulator3D( const FileReader & conf );

      /// Simulates a given time-length
      void simulate             ( real duration              );
      void simulateTimeStepCount();


      // Getter functions for the internally stored StaggeredGrid
      inline StaggeredGrid3D & grid();
      inline const StaggeredGrid3D & grid() const;

  private:
      std::string name;
      StaggeredGrid3D grid_;
      SORSolver3D solver_;
      real gamma, Re_1, tau;
      real GX, GY, GZ;
      int noOfTimeSteps, timeStepNumber, normFrequency;
      real time, dt;
      const real dtmax;
      real dx, dy, dz, dx_2, dy_2, dz_2;
      int outputInterval, imax, jmax, kmax;
      
      std::string boundary_condition_N, boundary_condition_S, boundary_condition_E,
                  boundary_condition_W, boundary_condition_U, boundary_condition_D;
      real boundary_velocity_N, boundary_velocity_S, boundary_velocity_E,
           boundary_velocity_W, boundary_velocity_U, boundary_velocity_D; 

      Array<real> &u = grid_.u();
      Array<real> &v = grid_.v();
      Array<real> &w = grid_.w();
      Array<real> &f = grid_.f();
      Array<real> &g = grid_.g();
      Array<real> &h = grid_.h();
      Array<real> &p = grid_.p();
      Array<real> &rhs = grid_.rhs();

      // Functions to compute derivatives
      inline real dU2_dx(int i, int j, int k);
      inline real dV2_dy(int i, int j, int k);
      inline real dW2_dz(int i, int j, int k);
      inline real dUV_dy(int i, int j, int k);
      inline real dUW_dz(int i, int j, int k);
      inline real dVU_dx(int i, int j, int k);
      inline real dVW_dz(int i, int j, int k);
      inline real dWU_dx(int i, int j, int k);
      inline real dWV_dy(int i, int j, int k);
      inline real d2U_dx2(int i, int j, int k);
      inline real d2U_dy2(int i, int j, int k);
      inline real d2U_dz2(int i, int j, int k);
      inline real d2V_dx2(int i, int j, int k);
      inline real d2V_dy2(int i, int j, int k);
      inline real d2V_dz2(int i, int j, int k);
      inline real d2W_dx2(int i, int j, int k);
      inline real d2W_dy2(int i, int j, int k);
      inline real d2W_dz2(int i, int j, int k);

	void refreshBoundaries();
	void normalizePressure();
      void computeFGH();
      void composeRHS();
      void determineNextDT( real const & limit);
      void updateVelocities();

};

inline StaggeredGrid3D & FluidSimulator3D::grid()
{
   return grid_;
}

inline const StaggeredGrid3D & FluidSimulator3D::grid() const
{
   return grid_;
}

#endif
