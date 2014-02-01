#ifndef __FLUID_SIMULATOR_H__
#define __FLUID_SIMULATOR_H__
#include "FileReader.hh"
#include "StaggeredGrid.hh"
#include "SORSolver.hh"
#include "RedBlackSORSolver.hh"
#include "CGSolver.hh"
#include "MultigridSolver.hh"
class FluidSimulator
{
  public:
      FluidSimulator( const FileReader & conf );

      /// Simulates a given time-length
      void simulate             ( real duration              );
      void simulateTimeStepCount();


      // Getter functions for the internally stored StaggeredGrid
      inline StaggeredGrid & grid();
      inline const StaggeredGrid & grid() const;
      ~FluidSimulator();

  private:
      std::string name;
      StaggeredGrid grid_;
      Solver* solver_;
      real gamma, Re_1, tau;
      real GX,GY;
      int noOfTimeSteps, timeStepNumber, normFrequency;
      real time, dt;
      const real dtmax;
      real dx, dy, dx_2, dy_2, inv4dx, inv4dy;
      int outputInterval, imax, jmax;
      
      std::string boundary_condition_N, boundary_condition_S, boundary_condition_E, boundary_condition_W;
      real boundary_velocity_N, boundary_velocity_S, boundary_velocity_E, boundary_velocity_W; 

      Array<real> &u = grid_.u();
      Array<real> &v = grid_.v();
      Array<real> &f = grid_.f();
      Array<real> &g = grid_.g();
      Array<real> &p = grid_.p();
      Array<real> &rhs = grid_.rhs();

      // Functions to compute derivatives
      inline real dU2_dx(int i, int j);
      inline real dV2_dy(int i, int j);
      inline real dUV_dy(int i, int j);
      inline real dVU_dx(int i, int j);
      inline real d2U_dx2(int i, int j);
      inline real d2U_dy2(int i, int j);
      inline real d2V_dx2(int i, int j);
      inline real d2V_dy2(int i, int j);

	void refreshBoundaries();
	void normalizePressure();
      void computeFG();
      void composeRHS();
      void determineNextDT( real const & limit);
      void updateVelocities();

};

inline StaggeredGrid & FluidSimulator::grid()
{
   return grid_;
}

inline const StaggeredGrid & FluidSimulator::grid() const
{
   return grid_;
}

#endif
