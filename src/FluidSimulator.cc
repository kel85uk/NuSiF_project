#include "FluidSimulator.hh"
#include "VTKWriter.hh"
#include <cmath>

FluidSimulator::FluidSimulator( const FileReader & conf ) : 
					  name(conf.getStringParameter("name")), grid_(StaggeredGrid(conf)),
					  gamma(conf.getRealParameter("gamma")),
					  Re_1(1.0/conf.getRealParameter("Re")), tau(conf.getRealParameter("safetyfactor")),
					  GX(conf.getRealParameter("GX")), GY(conf.getRealParameter("GY")),
					  noOfTimeSteps(conf.getIntParameter("timesteps")),timeStepNumber(0),time(0.0),
					  dt(conf.getRealParameter("dt")), dtmax(conf.getRealParameter("dtmax")),
					  dx(grid_.dx()), dy(grid_.dy()), dx_2(1.0/(dx*dx)), dy_2(1.0/(dy*dy)),
					  inv4dx(0.25/dx), inv4dy(0.25/dy),
					  imax(conf.getIntParameter("imax")), jmax(conf.getIntParameter("jmax")),
					  boundary_condition_N("noslip"), boundary_condition_S("noslip"),
					  boundary_condition_E("noslip"), boundary_condition_W("noslip"),
					  boundary_velocity_N(0), boundary_velocity_S(0), boundary_velocity_E(0),
					  boundary_velocity_W(0)
{
   if (conf.getStringParameter("solver") == "SOR")
       solver_ = new SORSolver(conf);
   else if (conf.getStringParameter("solver") == "RedBlackSOR")
       solver_ = new RedBlackSORSolver(conf);
   else if (conf.getStringParameter("solver") == "CG")
       solver_ = new CGSolver(conf);
   else if (conf.getStringParameter("solver") == "Multigrid")
       solver_ = new MultigridSolver(conf);
   else 
       CHECK_MSG(0, "Invalid Solver. Valid solvers are\n1. SOR\n2. RedBlackSOR\n3. CG\n4. Multigrid");
   CHECK_MSG((Re_1 > 0), "Reynolds number (Re) must be greater than 0");
   CHECK_MSG((dt > 0), "dt must be greater than 0");
   CHECK_MSG((tau > 0 && tau <= 1), "safetyfactor must lie between 0 and 1");
   CHECK_MSG((dtmax > 0), "dtmax must be greater than 0");
   if (conf.find("normalizationfrequency"))
      normFrequency = conf.getIntParameter("normalizationfrequency");
   else {
      WARN("Normalization frequency not defined. Using default value 100");
      normFrequency = 100; }
      
   CHECK_MSG(conf.find("outputinterval"), "Output Interval not defined. Not writing to file");
   outputInterval = conf.getIntParameter("outputinterval");
   
   if (conf.find("boundary_condition_N")) {
      boundary_condition_N = conf.getStringParameter("boundary_condition_N");
      CHECK_MSG ( boundary_condition_N != "outflow" || !(conf.find("boundary_velocity_N")),
                  "OUTFLOW boundary condition cannot have velocity specified" );
      CHECK_MSG ( boundary_condition_N != "inflow" || (conf.find("boundary_velocity_N")),
                  "INFLOW boundary condition must have velocity specified" ); }
   else
      boundary_condition_N = "noslip"; 
   if (conf.find("boundary_velocity_N"))
      boundary_velocity_N = conf.getRealParameter("boundary_velocity_N");

   if (conf.find("boundary_condition_S")) {
      boundary_condition_S = conf.getStringParameter("boundary_condition_S");
      CHECK_MSG ( boundary_condition_S != "outflow" || !(conf.find("boundary_velocity_S")),
                  "OUTFLOW boundary condition cannot have velocity specified" );
      CHECK_MSG ( boundary_condition_S != "inflow" || (conf.find("boundary_velocity_S")),
                  "INFLOW boundary condition must have velocity specified" ); }
   else
      boundary_condition_S = "noslip"; 
   if (conf.find("boundary_velocity_S"))
      boundary_velocity_S = conf.getRealParameter("boundary_velocity_S");

   if (conf.find("boundary_condition_E")) {
      boundary_condition_E = conf.getStringParameter("boundary_condition_E");
      CHECK_MSG ( boundary_condition_E != "outflow" || !(conf.find("boundary_velocity_E")),
                  "OUTFLOW boundary condition cannot have velocity specified" );
      CHECK_MSG ( boundary_condition_E != "inflow" || (conf.find("boundary_velocity_E")),
                  "INFLOW boundary condition must have velocity specified" );
      }
   else
      boundary_condition_E = "noslip"; 
   if (conf.find("boundary_velocity_E"))
      boundary_velocity_E = conf.getRealParameter("boundary_velocity_E");
      
   if (conf.find("boundary_condition_W")) {
      boundary_condition_W = conf.getStringParameter("boundary_condition_W");
      CHECK_MSG ( boundary_condition_W != "outflow" || !(conf.find("boundary_velocity_W")),
                  "OUTFLOW boundary condition cannot have velocity specified" );
      CHECK_MSG ( boundary_condition_W != "inflow" || (conf.find("boundary_velocity_W")),
                  "INFLOW boundary condition must have velocity specified" );
      }
   else
      boundary_condition_W = "noslip"; 
   if (conf.find("boundary_velocity_W"))
      boundary_velocity_W = conf.getRealParameter("boundary_velocity_W");
   
}

FluidSimulator::~FluidSimulator()
{
   delete [] solver_;
}

inline real FluidSimulator::dU2_dx(int i, int j)
{
   return (inv4dx * (  pow(grid_.u(i,j,CENTER)+grid_.u(i,j,EAST),2) - pow(grid_.u(i,j,WEST)+grid_.u(i,j,CENTER),2) 
                             + gamma * (  std::abs(grid_.u(i,j,CENTER)+grid_.u(i,j,EAST))
                                          * (grid_.u(i,j,CENTER)-grid_.u(i,j,EAST))
                                        - std::abs(grid_.u(i,j,WEST)+grid_.u(i,j,CENTER))
                                          * (grid_.u(i,j,WEST)-grid_.u(i,j,CENTER)) ) ) );
}
inline real FluidSimulator::dV2_dy(int i, int j)
{
   return (inv4dy * (  pow(grid_.v(i,j,CENTER)+grid_.v(i,j,NORTH),2) - pow(grid_.v(i,j,SOUTH)+grid_.v(i,j,CENTER),2) 
                             + gamma * (  std::abs(grid_.v(i,j,CENTER)+grid_.v(i,j,NORTH))
                                          * (grid_.v(i,j,CENTER)-grid_.v(i,j,NORTH)) 
                                        - std::abs(grid_.v(i,j,SOUTH)+grid_.v(i,j,CENTER))
                                          * (grid_.v(i,j,SOUTH)-grid_.v(i,j,CENTER)) ) ) );
}
inline real FluidSimulator::dUV_dy(int i, int j)
{
   return (inv4dy * (  (grid_.v(i,j,CENTER)+grid_.v(i,j,EAST))*(grid_.u(i,j,CENTER)+grid_.u(i,j,NORTH))
                      - (grid_.v(i,j,SOUTH)+grid_.v(i+1,j,SOUTH))*(grid_.u(i,j,SOUTH)+grid_.u(i,j,CENTER))
                      + gamma * (  std::abs(grid_.v(i,j,CENTER)+grid_.v(i,j,EAST)) * (grid_.u(i,j,CENTER)-grid_.u(i,j,NORTH))  
                                 - std::abs(grid_.v(i,j,SOUTH)+grid_.v(i+1,j,SOUTH)) * (grid_.u(i,j,SOUTH)-grid_.u(i,j,CENTER)) ) ) );
}
inline real FluidSimulator::dVU_dx(int i, int j)
{
   return (inv4dx * (  (grid_.u(i,j,CENTER)+grid_.u(i,j,NORTH))*(grid_.v(i,j,CENTER)+grid_.v(i,j,EAST))
                      - (grid_.u(i,j,WEST)+grid_.u(i,j+1,WEST))*(grid_.v(i,j,WEST)+grid_.v(i,j,CENTER)) 
                      + gamma * (  std::abs(grid_.u(i,j,CENTER)+grid_.u(i,j,NORTH)) * (grid_.v(i,j,CENTER)-grid_.v(i,j,EAST)) 
                                 - std::abs(grid_.u(i,j,WEST)+grid_.u(i,j+1,WEST)) * (grid_.v(i,j,WEST)-grid_.v(i,j,CENTER)) ) ) );
}
inline real FluidSimulator::d2U_dx2(int i, int j)
{
   return (dx_2 * ( grid_.u(i,j,EAST)-2.0*grid_.u(i,j,CENTER)+grid_.u(i,j,WEST) ) );
}
inline real FluidSimulator::d2U_dy2(int i, int j)
{
   return (dy_2 * ( grid_.u(i,j,NORTH)-2.0*grid_.u(i,j,CENTER)+grid_.u(i,j,SOUTH) ) );
}
inline real FluidSimulator::d2V_dx2(int i, int j)
{
   return (dx_2 * ( grid_.v(i,j,EAST)-2.0*grid_.v(i,j,CENTER)+grid_.v(i,j,WEST) ) );
}
inline real FluidSimulator::d2V_dy2(int i, int j)
{
   return (dy_2 * ( grid_.v(i,j,NORTH)-2.0*grid_.v(i,j,CENTER)+grid_.v(i,j,SOUTH) ) );
}

void FluidSimulator::determineNextDT( real const & limit)
{
   real stepSize = tau*std::min( (0.5/Re_1) / (dx_2+dy_2),
                                  std::min( dx/std::abs(u.maximum()), dy/std::abs(v.maximum()) ) );
   dt = stepSize<limit ? stepSize : limit;
}
void FluidSimulator::refreshBoundaries()
{
   if (boundary_condition_N == "noslip")
       for (int i=1; i<=imax; i++) {
           u(i,jmax + 1) = 2.0 * boundary_velocity_N - u(i,jmax);
           v(i,jmax) = 0.0; }
   else if (boundary_condition_N == "inflow")
       for (int i=1; i<=imax; i++) {
           u(i,jmax + 1) = -u(i,jmax);
           v(i,jmax) = - boundary_velocity_N; }
   else if (boundary_condition_N == "outflow")
       for (int i=1; i<=imax; i++) {
           u(i,jmax + 1) = u(i,jmax);
           v(i,jmax) = v(i,jmax-1); }
           
   if (boundary_condition_S == "noslip")
       for (int i=1; i<=imax; i++) {
           u(i,0) = 2.0 * boundary_velocity_S -u(i,1);
           v(i,0) = 0.0; }
   else if (boundary_condition_S == "inflow")
       for (int i=1; i<=imax; i++) {
           u(i,0) = -u(i,1);
           v(i,0) = boundary_velocity_S; }
   else if (boundary_condition_S == "outflow")
       for (int i=1; i<=imax; i++) {
           u(i,0) = u(i,1);
           v(i,0) = v(i,1); }
          
   if (boundary_condition_E == "noslip")
       for (int j=1; j<=jmax; j++) {
           u(imax, j) = 0.0;
           v(imax + 1, j) = 2.0 * boundary_velocity_E -v(imax,j); }
   else if (boundary_condition_E == "inflow")
       for (int j=1; j<=jmax; j++) {
           u(imax, j) = -boundary_velocity_E;
           v(imax + 1, j) = -v(imax,j); }
   else if (boundary_condition_E == "outflow")
       for (int j=1; j<=jmax; j++) {
           u(imax, j) = u(imax-1, j);
           v(imax + 1, j) = v(imax,j); }

   if (boundary_condition_W == "noslip")
       for (int j=1; j<=jmax; j++) {
           u(0,j) = 0.0;
           v(0,j) = 2.0 * boundary_velocity_W -v(1,j); }
   else if (boundary_condition_W == "inflow")
       for (int j=1; j<=jmax; j++) {
           u(0,j) = boundary_velocity_W;
           v(0,j) = -v(1,j); }
   else if (boundary_condition_W == "outflow")
       for (int j=1; j<=jmax; j++) {
           u(0,j) = u(1,j);
           v(0,j) = v(1,j); }     
}

void FluidSimulator::normalizePressure()
{
   real Pmean = 0.0;
   for (int i=1; i<=imax; i++) 
       for (int j=1; j<=jmax; j++)
           if (grid_.isFluid(i,j))
               Pmean += p(i,j);
   Pmean /= grid_.getNumFluid();
   for (int i=1; i<=imax; i++) 
       for (int j=1; j<=jmax; j++)
           if (grid_.isFluid(i,j))
               p(i,j) -= Pmean;
}

void FluidSimulator::computeFG()
{
   for (int i=1; i<imax; i++) 
       for (int j=1; j<=jmax; j++)
           f(i,j) = grid_.u(i,j,CENTER) + dt *  ( Re_1*(d2U_dx2(i,j)+d2U_dy2(i,j)) - dU2_dx(i,j) - dUV_dy(i,j) + GX);

   for (int i=1; i<=imax; i++) 
       for (int j=1; j<jmax; j++)
           g(i,j) = grid_.v(i,j,CENTER) + dt *  ( Re_1*(d2V_dx2(i,j)+d2V_dy2(i,j)) - dV2_dy(i,j) - dVU_dx(i,j) + GY);

   for (int i=0; i<=imax; i++){
       g(i,0) = v(i,0);
       g(i,jmax) = v(i,jmax); }        

   for (int j=0; j<=jmax; j++) {
       f(0,j) = u(0,j);
       f(imax,j) = u(imax,j); }           
}

void FluidSimulator::composeRHS()
{
   for (int i=1; i<=imax; i++) 
       for (int j=1; j<=jmax; j++)
           if (grid_.isFluid(i,j))
              rhs(i,j) = (1.0/dt) * ( (grid_.f(i,j,CENTER) - grid_.f(i,j,WEST))/dx + (grid_.g(i,j,CENTER) - grid_.g(i,j,SOUTH))/dy );
}

void FluidSimulator::updateVelocities()
{
   for (int i=1; i<imax; i++) 
       for (int j=1; j<=jmax; j++)
           if (grid_.isFluid(i,j) && grid_.isFluid(i+1,j))
               u(i,j) = f(i,j) - dt/dx *  ( p(i+1,j) - p(i,j) );
   for (int i=1; i<=imax; i++) 
       for (int j=1; j<jmax; j++)
           if (grid_.isFluid(i,j) && grid_.isFluid(i,j+1))
               v(i,j) = g(i,j) - dt/dy *  ( p(i,j+1) - p(i,j) );
}

void FluidSimulator::simulateTimeStepCount()
{
   VTKWriter vtkWriter(grid_, name, true, true);
   while (timeStepNumber<=noOfTimeSteps)
   {
      if (tau>0)
         determineNextDT( dtmax );
      refreshBoundaries();
      if (timeStepNumber%outputInterval==0)
         vtkWriter.write();
      if (timeStepNumber%normFrequency == 0)
         normalizePressure();
      computeFG();
      composeRHS();
      std::cout<< "\nTime step    = "<<timeStepNumber << "\tTime     = "<< time<<"\n";
      solver_->solve(grid_);   
      updateVelocities();
      timeStepNumber++;
      time += dt;
   }
}

void FluidSimulator::simulate( real duration )
{
   VTKWriter vtkWriter(grid_, name, true, true);
   while (time<=duration)
   {
      if (tau>0)
         determineNextDT( dtmax );
      refreshBoundaries();
      if (timeStepNumber%outputInterval==0)
         vtkWriter.write();
      if (timeStepNumber%normFrequency == 0)
         normalizePressure();
      computeFG();
      composeRHS();
      std::cout<< "\nTime step    = "<<timeStepNumber << "\tTime     = "<< time<<"\n";
      solver_->solve(grid_);   
      updateVelocities();
      timeStepNumber++;
      time += dt;
   }
}
