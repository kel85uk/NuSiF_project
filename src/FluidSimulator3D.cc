#include "FluidSimulator3D.hh"
#include "VTKWriter3D.hh"
#include <cmath>

FluidSimulator3D::FluidSimulator3D( const FileReader & conf ) : 
					  name(conf.getStringParameter("name")), grid_(StaggeredGrid3D(conf)),
					  gamma(conf.getRealParameter("gamma")),
					  Re_1(1.0/conf.getRealParameter("Re")), tau(conf.getRealParameter("safetyfactor")),
					  GX(conf.getRealParameter("GX")), GY(conf.getRealParameter("GY")), GZ(conf.getRealParameter("GZ")),
					  noOfTimeSteps(conf.getIntParameter("timesteps")),timeStepNumber(0),time(0.0),
					  dt(conf.getRealParameter("dt")), dtmax(conf.getRealParameter("dtmax")),
					  dx(grid_.dx()), dy(grid_.dy()), dz(grid_.dz()),
					  dx_2(1.0/(dx*dx)), dy_2(1.0/(dy*dy)), dz_2(1.0/(dz*dz)),
					  imax(conf.getIntParameter("imax")), jmax(conf.getIntParameter("jmax")),
					  kmax(conf.getIntParameter("kmax")),
					  boundary_condition_N("noslip"), boundary_condition_S("noslip"),
					  boundary_condition_E("noslip"), boundary_condition_W("noslip"),
					  boundary_condition_U("noslip"), boundary_condition_D("noslip"),
					  boundary_velocity_N1(0), boundary_velocity_N2(0), boundary_velocity_N3(0),
					  boundary_velocity_S1(0), boundary_velocity_S2(0), boundary_velocity_S3(0),
					  boundary_velocity_E1(0), boundary_velocity_E2(0), boundary_velocity_E3(0),
					  boundary_velocity_W1(0), boundary_velocity_W2(0), boundary_velocity_W3(0),
					  boundary_velocity_U1(0), boundary_velocity_U2(0), boundary_velocity_U3(0),
					  boundary_velocity_D1(0), boundary_velocity_D2(0), boundary_velocity_D3(0)
{
   if (conf.getStringParameter("solver") == "SOR")
       solver_ = new SORSolver3D(conf);
   else if (conf.getStringParameter("solver") == "RedBlackSOR")
       solver_ = new RedBlackSORSolver3D(conf);
   else if (conf.getStringParameter("solver") == "CG")
       solver_ = new CGSolver3D(conf);
/*   else if (conf.getStringParameter("solver") == "Multigrid")
       solver_ = new MultigridSolver3D(conf);
*/   else 
       CHECK_MSG(0, "Invalid Solver. Valid solvers are/n1. SOR/n2. RedBlackSOR\n3. CG\n4. Multigrid");
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

      CHECK_MSG ((boundary_condition_N != "outflow")
                  || !(conf.find("boundary_velocity_N1") || conf.find("boundary_velocity_N2") || conf.find("boundary_velocity_N3")),
                 "outflow boundary condition cannot have velocity specified");  
      CHECK_MSG ( boundary_condition_N != "inflow" || !(conf.find("boundary_velocity_N1") || conf.find("boundary_velocity_N3")),
                  "inflow boundary condition must only have normal velocity specified" );
      CHECK_MSG ( boundary_condition_N != "inflow" || (conf.find("boundary_velocity_N2")),
                  "inflow boundary condition must have normal velocity specified" ); }
   else
      boundary_condition_N = "noslip"; 
   if (conf.find("boundary_velocity_N1"))
      boundary_velocity_N1 = conf.getRealParameter("boundary_velocity_N1");
   if (conf.find("boundary_velocity_N3"))
      boundary_velocity_N3 = conf.getRealParameter("boundary_velocity_N3");

   if (conf.find("boundary_condition_S")) {
      boundary_condition_S = conf.getStringParameter("boundary_condition_S");
      CHECK_MSG ((boundary_condition_S != "outflow")
                  || !(conf.find("boundary_velocity_S1") || conf.find("boundary_velocity_S2") || conf.find("boundary_velocity_S3")),
                 "outflow boundary condition cannot have velocity specified");  
      CHECK_MSG ( boundary_condition_S != "inflow" || !(conf.find("boundary_velocity_S1") || conf.find("boundary_velocity_S3")),
                  "inflow boundary condition must only have normal velocity specified" );
      CHECK_MSG ( boundary_condition_S != "inflow" || (conf.find("boundary_velocity_S2")),
                  "inflow boundary condition must have normal velocity specified" ); }
   else
      boundary_condition_S = "noslip"; 
   if (conf.find("boundary_velocity_S1"))
      boundary_velocity_S1 = conf.getRealParameter("boundary_velocity_S1");
   if (conf.find("boundary_velocity_S3"))
      boundary_velocity_S3 = conf.getRealParameter("boundary_velocity_S3");

   if (conf.find("boundary_condition_E")) {
      boundary_condition_E = conf.getStringParameter("boundary_condition_E");
      CHECK_MSG ((boundary_condition_E != "outflow")
                  || !(conf.find("boundary_velocity_E1") || conf.find("boundary_velocity_E2") || conf.find("boundary_velocity_E3")),
                 "outflow boundary condition cannot have velocity specified");  
      CHECK_MSG ( boundary_condition_E != "inflow" || !(conf.find("boundary_velocity_E2") || conf.find("boundary_velocity_E3")),
                  "inflow boundary condition must only have normal velocity specified" );
      CHECK_MSG ( boundary_condition_E != "inflow" || (conf.find("boundary_velocity_E1")),
                  "inflow boundary condition must have normal velocity specified" ); }
   else
      boundary_condition_E = "noslip"; 
   if (conf.find("boundary_velocity_E2"))
      boundary_velocity_E2 = conf.getRealParameter("boundary_velocity_E2");
   if (conf.find("boundary_velocity_E3"))
      boundary_velocity_E3 = conf.getRealParameter("boundary_velocity_E3");
      
   if (conf.find("boundary_condition_W")) {
      boundary_condition_W = conf.getStringParameter("boundary_condition_W");
      CHECK_MSG ((boundary_condition_W != "outflow")
                  || !(conf.find("boundary_velocity_W1") || conf.find("boundary_velocity_W2") || conf.find("boundary_velocity_W3")),
                 "outflow boundary condition cannot have velocity specified");  
      CHECK_MSG ( boundary_condition_W != "inflow" || !(conf.find("boundary_velocity_W2") || conf.find("boundary_velocity_W3")),
                  "inflow boundary condition must only have normal velocity specified" );
      CHECK_MSG ( boundary_condition_W != "inflow" || (conf.find("boundary_velocity_W1")),
                  "inflow boundary condition must have normal velocity specified" ); }
   else
      boundary_condition_W = "noslip"; 
   if (conf.find("boundary_velocity_W2"))
      boundary_velocity_W2 = conf.getRealParameter("boundary_velocity_W2");
   if (conf.find("boundary_velocity_W3"))
      boundary_velocity_W3 = conf.getRealParameter("boundary_velocity_W3");
   
   if (conf.find("boundary_condition_U")) {
      boundary_condition_U = conf.getStringParameter("boundary_condition_U");
      CHECK_MSG ((boundary_condition_U != "outflow")
                  || !(conf.find("boundary_velocity_U1") || conf.find("boundary_velocity_U2") || conf.find("boundary_velocity_U3")),
                 "outflow boundary condition cannot have velocity specified");  
      CHECK_MSG ( boundary_condition_U != "inflow" || !(conf.find("boundary_velocity_U1") || conf.find("boundary_velocity_U2")),
                  "inflow boundary condition must only have normal velocity specified" );
      CHECK_MSG ( boundary_condition_U != "inflow" || (conf.find("boundary_velocity_U3")),
                  "inflow boundary condition must have normal velocity specified" ); }
   else
      boundary_condition_U = "noslip"; 
   if (conf.find("boundary_velocity_U1"))
      boundary_velocity_U1 = conf.getRealParameter("boundary_velocity_U1");
   if (conf.find("boundary_velocity_U2"))
      boundary_velocity_U2 = conf.getRealParameter("boundary_velocity_U2");
   
   if (conf.find("boundary_condition_D")) {
      boundary_condition_D = conf.getStringParameter("boundary_condition_D");
      CHECK_MSG ((boundary_condition_D != "outflow")
                  || !(conf.find("boundary_velocity_D1") || conf.find("boundary_velocity_D2") || conf.find("boundary_velocity_D3")),
                 "outflow boundary condition cannot have velocity specified");  
      CHECK_MSG ( boundary_condition_D != "inflow" || !(conf.find("boundary_velocity_D1") || conf.find("boundary_velocity_D2")),
                  "inflow boundary condition must only have normal velocity specified" );
      CHECK_MSG ( boundary_condition_D != "inflow" || (conf.find("boundary_velocity_D3")),
                  "inflow boundary condition must have normal velocity specified" ); }
   else
      boundary_condition_D = "noslip"; 
   if (conf.find("boundary_velocity_D1"))
      boundary_velocity_D1 = conf.getRealParameter("boundary_velocity_D1");
   if (conf.find("boundary_velocity_D2"))
      boundary_velocity_D2 = conf.getRealParameter("boundary_velocity_D2");
   
}

FluidSimulator3D::~FluidSimulator3D()
{
   delete [] solver_;
}

inline real FluidSimulator3D::dU2_dx(int i, int j, int k)
{
return (0.25/dx*(pow(grid_.u(i,j,k,CENTER)+grid_.u(i,j,k,EAST),2) - pow(grid_.u(i,j,k,WEST)+grid_.u(i,j,k,CENTER),2) 
                    +gamma*( std::abs(grid_.u(i,j,k,CENTER)+grid_.u(i,j,k,EAST)) * (grid_.u(i,j,k,CENTER)-grid_.u(i,j,k,EAST)) 
                            -std::abs(grid_.u(i,j,k,WEST)+grid_.u(i,j,k,CENTER)) * (grid_.u(i,j,k,WEST)-grid_.u(i,j,k,CENTER)))));
}

inline real FluidSimulator3D::dV2_dy(int i, int j, int k)
{
return (0.25/dy*(pow(grid_.v(i,j,k,CENTER)+grid_.v(i,j,k,NORTH),2) - pow(grid_.v(i,j,k,SOUTH)+grid_.v(i,j,k,CENTER),2) 
                    +gamma*( std::abs(grid_.v(i,j,k,CENTER)+grid_.v(i,j,k,NORTH)) * (grid_.v(i,j,k,CENTER)-grid_.v(i,j,k,NORTH)) 
                            -std::abs(grid_.v(i,j,k,SOUTH)+grid_.v(i,j,k,CENTER)) * (grid_.v(i,j,k,SOUTH)-grid_.v(i,j,k,CENTER)))));
}

inline real FluidSimulator3D::dW2_dz(int i, int j, int k)
{
return (0.25/dy*(pow(grid_.w(i,j,k,CENTER)+grid_.w(i,j,k,UP),2) - pow(grid_.w(i,j,k,CENTER)+grid_.w(i,j,k,DOWN),2) 
                     + gamma*( std::abs(grid_.w(i,j,k,CENTER)+grid_.w(i,j,k,UP)) * (grid_.w(i,j,k,CENTER)-grid_.w(i,j,k,UP)) 
                              -std::abs(grid_.w(i,j,k,DOWN)+grid_.w(i,j,k,CENTER)) * (grid_.w(i,j,k,CENTER)-grid_.w(i,j,k,DOWN)))));
}

inline real FluidSimulator3D::dUV_dy(int i, int j, int k)
{
return (0.25/dy* ( (grid_.v(i,j,k,CENTER)+grid_.v(i,j,k,EAST))*(grid_.u(i,j,k,CENTER)+grid_.u(i,j,k,NORTH))
                  -(grid_.v(i,j,k,SOUTH)+grid_.v(i+1,j,k,SOUTH))*(grid_.u(i,j,k,SOUTH)+grid_.u(i,j,k,CENTER))
                  +gamma * (  std::abs(grid_.v(i,j,k,CENTER)+grid_.v(i,j,k,EAST)) * (grid_.u(i,j,k,CENTER)-grid_.u(i,j,k,NORTH))  
                            - std::abs(grid_.v(i,j,k,SOUTH)+grid_.v(i+1,j,k,SOUTH)) * (grid_.u(i,j,k,SOUTH)-grid_.u(i,j,k,CENTER)))));
}
inline real FluidSimulator3D::dUW_dz(int i, int j, int k)
{
return (0.25/dz* ( (grid_.w(i,j,k,CENTER)+grid_.w(i,j,k,EAST))*(grid_.u(i,j,k,CENTER)+grid_.u(i,j,k,UP))
                  -(grid_.w(i,j,k,DOWN)+grid_.w(i+1,j,k,DOWN))*(grid_.u(i,j,k,DOWN)+grid_.u(i,j,k,CENTER))
                  +gamma * (  std::abs(grid_.w(i,j,k,CENTER)+grid_.w(i,j,k,EAST)) * (grid_.u(i,j,k,CENTER)-grid_.u(i,j,k,UP))  
                            - std::abs(grid_.w(i,j,k,DOWN)+grid_.w(i+1,j,k,DOWN)) * (grid_.u(i,j,k,DOWN)-grid_.u(i,j,k,CENTER)))));
}
inline real FluidSimulator3D::dVU_dx(int i, int j, int k)
{
return (0.25/dx * (  (grid_.u(i,j,k,CENTER)+grid_.u(i,j,k,NORTH))*(grid_.v(i,j,k,CENTER)+grid_.v(i,j,k,EAST))
                   - (grid_.u(i,j,k,WEST)+grid_.u(i,j+1,k,WEST))*(grid_.v(i,j,k,WEST)+grid_.v(i,j,k,CENTER)) 
                   + gamma * (  std::abs(grid_.u(i,j,k,CENTER)+grid_.u(i,j,k,NORTH)) * (grid_.v(i,j,k,CENTER)-grid_.v(i,j,k,EAST)) 
                              - std::abs(grid_.u(i,j,k,WEST)+grid_.u(i,j+1,k,WEST)) * (grid_.v(i,j,k,WEST)-grid_.v(i,j,k,CENTER)))));
}
inline real FluidSimulator3D::dVW_dz(int i, int j, int k)
{
return (0.25/dz * (  (grid_.w(i,j,k,CENTER)+grid_.w(i,j,k,NORTH))*(grid_.v(i,j,k,CENTER)+grid_.v(i,j,k,UP))
                   - (grid_.w(i,j,k,DOWN)+grid_.w(i,j+1,k,DOWN))*(grid_.v(i,j,k,DOWN)+grid_.v(i,j,k,CENTER)) 
                   + gamma * (  std::abs(grid_.w(i,j,k,CENTER)+grid_.w(i,j,k,NORTH)) * (grid_.v(i,j,k,CENTER)-grid_.v(i,j,k,UP)) 
                              - std::abs(grid_.w(i,j,k,DOWN)+grid_.w(i,j+1,k,DOWN)) * (grid_.v(i,j,k,DOWN)-grid_.v(i,j,k,CENTER)))));
}
inline real FluidSimulator3D::dWU_dx(int i, int j, int k)
{
return (0.25/dx * (  (grid_.u(i,j,k,CENTER)+grid_.u(i,j,k,UP))*(grid_.w(i,j,k,CENTER)+grid_.w(i,j,k,EAST))
                   - (grid_.u(i,j,k,WEST)+grid_.u(i,j,k+1,WEST))*(grid_.w(i,j,k,WEST)+grid_.w(i,j,k,CENTER)) 
                   + gamma * (  std::abs(grid_.u(i,j,k,CENTER)+grid_.u(i,j,k,UP)) * (grid_.w(i,j,k,CENTER)-grid_.w(i,j,k,EAST)) 
                              - std::abs(grid_.u(i,j,k,WEST)+grid_.u(i,j,k+1,WEST)) * (grid_.w(i,j,k,WEST)-grid_.w(i,j,k,CENTER)))));
}
inline real FluidSimulator3D::dWV_dy(int i, int j, int k)
{
return (0.25/dy* ( (grid_.v(i,j,k,CENTER)+grid_.v(i,j,k,UP))*(grid_.w(i,j,k,CENTER)+grid_.w(i,j,k,NORTH))
                  -(grid_.v(i,j,k,SOUTH)+grid_.v(i,j,k+1,SOUTH))*(grid_.w(i,j,k,SOUTH)+grid_.w(i,j,k,CENTER))
                  +gamma * (  std::abs(grid_.v(i,j,k,CENTER)+grid_.v(i,j,k,UP)) * (grid_.w(i,j,k,CENTER)-grid_.w(i,j,k,NORTH))  
                            - std::abs(grid_.v(i,j,k,SOUTH)+grid_.v(i,j,k+1,SOUTH)) * (grid_.w(i,j,k,SOUTH)-grid_.w(i,j,k,CENTER)))));
}

inline real FluidSimulator3D::d2U_dx2(int i, int j, int k)
{
   return (dx_2 * ( grid_.u(i,j,k,EAST)-2.0*grid_.u(i,j,k,CENTER)+grid_.u(i,j,k,WEST) ) );
}
inline real FluidSimulator3D::d2U_dy2(int i, int j, int k)
{
   return (dy_2 * ( grid_.u(i,j,k,NORTH)-2.0*grid_.u(i,j,k,CENTER)+grid_.u(i,j,k,SOUTH) ) );
}
inline real FluidSimulator3D::d2U_dz2(int i, int j, int k)
{
   return (dz_2 * ( grid_.u(i,j,k,UP)-2.0*grid_.u(i,j,k,CENTER)+grid_.u(i,j,k,DOWN) ) );
}
inline real FluidSimulator3D::d2V_dx2(int i, int j, int k)
{
   return (dx_2 * ( grid_.v(i,j,k,EAST)-2.0*grid_.v(i,j,k,CENTER)+grid_.v(i,j,k,WEST) ) );
}
inline real FluidSimulator3D::d2V_dy2(int i, int j, int k)
{
   return (dy_2 * ( grid_.v(i,j,k,NORTH)-2.0*grid_.v(i,j,k,CENTER)+grid_.v(i,j,k,SOUTH) ) );
}
inline real FluidSimulator3D::d2V_dz2(int i, int j, int k)
{
   return (dz_2 * ( grid_.v(i,j,k,UP)-2.0*grid_.v(i,j,k,CENTER)+grid_.v(i,j,k,DOWN) ) );
}
inline real FluidSimulator3D::d2W_dx2(int i, int j, int k)
{
   return (dx_2 * ( grid_.w(i,j,k,EAST)-2.0*grid_.w(i,j,k,CENTER)+grid_.w(i,j,k,WEST) ) );
}
inline real FluidSimulator3D::d2W_dy2(int i, int j, int k)
{
   return (dy_2 * ( grid_.w(i,j,k,NORTH)-2.0*grid_.w(i,j,k,CENTER)+grid_.w(i,j,k,SOUTH) ) );
}
inline real FluidSimulator3D::d2W_dz2(int i, int j, int k)
{
   return (dz_2 * ( grid_.w(i,j,k,UP)-2.0*grid_.w(i,j,k,CENTER)+grid_.w(i,j,k,DOWN) ) );
}

void FluidSimulator3D::determineNextDT( real const & limit)
{
   real stepSize = tau*std::min( (0.5/Re_1) / (dx_2+dy_2+dz_2),
                                  std::min( dx/std::abs(u.maximum()),
                                            std::min( dy/std::abs(v.maximum()), dz/std::abs(w.maximum()) ) ) );
   dt = stepSize<limit ? stepSize : limit;
}
void FluidSimulator3D::refreshBoundaries()
{
   if (boundary_condition_N == "noslip")
       for (int i=1; i<=imax; i++)
           for (int k=1; k<=kmax; k++) {
               u(i,jmax+1,k) = 2.0 * boundary_velocity_N1 - u(i,jmax,k);
               v(i,jmax,k) = 0.0;
               w(i,jmax+1,k) =  2.0 * boundary_velocity_N3 - w(i,jmax,k); }
   else if (boundary_condition_N == "inflow")
       for (int i=1; i<=imax; i++)
           for (int k=1; k<=kmax; k++) {
               u(i,jmax+1,k) = -u(i,jmax,k);
               v(i,jmax,k) = - boundary_velocity_N2;
               w(i,jmax+1,k) = -w(i,jmax,k); }
   else if (boundary_condition_N == "outflow")
       for (int i=1; i<=imax; i++)
           for (int k=1; k<=kmax; k++) {
               u(i,jmax+1,k) = u(i,jmax,k);
               v(i,jmax,k) = v(i,jmax-1,k);
               w(i,jmax+1,k) = w(i,jmax,k); }
           
   if (boundary_condition_S == "noslip")
       for (int i=1; i<=imax; i++)
           for (int k=1; k<=kmax; k++) {
               u(i,0,k) = 2.0 * boundary_velocity_S1 -u(i,1,k);
               v(i,0,k) = 0.0;
               w(i,0,k) = 2.0 * boundary_velocity_S3 -w(i,1,k); }
   else if (boundary_condition_S == "inflow")
       for (int i=1; i<=imax; i++)
           for (int k=1; k<=kmax; k++) {
               u(i,0,k) = -u(i,1,k);
               v(i,0,k) = boundary_velocity_S2;
               w(i,0,k) = -w(i,1,k); }
   else if (boundary_condition_S == "outflow")
       for (int i=1; i<=imax; i++)
           for (int k=1; k<=kmax; k++) {
               u(i,0,k) = u(i,1,k);
               v(i,0,k) = v(i,1,k);
               w(i,0,k) = w(i,1,k); }
           
   if (boundary_condition_E == "noslip")
       for (int j=1; j<=jmax; j++)
           for (int k=1; k<=kmax; k++) {
               u(imax,j,k) = 0.0;
               v(imax+1,j,k) = 2.0 * boundary_velocity_E2 -v(imax,j,k);
               w(imax+1,j,k) = 2.0 * boundary_velocity_E3 -w(imax,j,k); }
   else if (boundary_condition_E == "inflow")
       for (int j=1; j<=jmax; j++)
           for (int k=1; k<=kmax; k++) {
               u(imax,j,k) = -boundary_velocity_E1;
               v(imax+1,j,k) = -v(imax,j,k);
               w(imax+1,j,k) = -w(imax,j,k); }
   else if (boundary_condition_E == "outflow")
       for (int j=1; j<=jmax; j++)
           for (int k=1; k<=kmax; k++) {
               u(imax,j,k) = u(imax-1,j,k);
               v(imax+1,j,k) = v(imax,j,k);
               w(imax+1,j,k) = w(imax,j,k); }
  
   if (boundary_condition_W == "noslip")
       for (int j=1; j<=jmax; j++)
           for (int k=1; k<=kmax; k++) {
               u(0,j,k) = 0.0;
               v(0,j,k) = 2.0 * boundary_velocity_W2 -v(1,j,k);
               w(0,j,k) = 2.0 * boundary_velocity_W3 -w(1,j,k); }
   else if (boundary_condition_W == "inflow")
       for (int j=1; j<=jmax; j++)
           for (int k=1; k<=kmax; k++) {
               u(0,j,k) = boundary_velocity_W1;
               v(0,j,k) = -v(1,j,k);
               w(0,j,k) = -w(1,j,k); }
   else if (boundary_condition_W == "outflow")
       for (int j=1; j<=jmax; j++)
           for (int k=1; k<=kmax; k++) {
               u(0,j,k) = u(1,j,k);
               v(0,j,k) = v(1,j,k);    
               w(0,j,k) = w(1,j,k); }     

   if (boundary_condition_U == "noslip")
       for (int i=1; i<=imax; i++)
           for (int j=1; j<=jmax; j++) {
               u(i,j,kmax+1) = 2.0 * boundary_velocity_U1 -u(i,j,kmax);
               v(i,j,kmax+1) = 2.0 * boundary_velocity_U2 -v(i,j,kmax);
               w(i,j,kmax) = 0.0; }
   else if (boundary_condition_U == "inflow")
       for (int i=1; i<=imax; i++)
           for (int j=1; j<=jmax; j++) {
               u(i,j,kmax+1) = -u(i,j,kmax);
               v(i,j,kmax+1) = -v(i,j,kmax);
               w(i,j,kmax) = -boundary_velocity_U3; }
   else if (boundary_condition_U == "outflow")
       for (int i=1; i<=imax; i++)
           for (int j=1; j<=jmax; j++) {
               u(i,j,kmax+1) = u(i,j,kmax);
               v(i,j,kmax+1) = v(i,j,kmax);
               w(i,j,kmax) = w(i,j,kmax-1); }
  
   if (boundary_condition_D == "noslip")
       for (int i=1; i<=imax; i++)
           for (int j=1; j<=jmax; j++) {
               u(i,j,0) = 2.0 * boundary_velocity_D1 -u(i,j,1);
               v(i,j,0) = 2.0 * boundary_velocity_D2 -v(i,j,1);
               w(i,j,0) = 0.0; }
   else if (boundary_condition_D == "inflow")
       for (int i=1; i<=imax; i++)
           for (int j=1; j<=jmax; j++) {
               u(i,j,0) = -u(i,j,1);
               v(i,j,0) = -v(i,j,1);
               w(i,j,0) = boundary_velocity_D3; }
   else if (boundary_condition_D == "outflow")
       for (int i=1; i<=imax; i++)
           for (int j=1; j<=jmax; j++) {
               u(i,j,0) = u(i,j,0);
               v(i,j,0) = v(i,j,0);    
               w(i,j,0) = w(i,j,0); }     
}

void FluidSimulator3D::normalizePressure()
{
   real Pmean = 0.0;
   for (int i=1; i<=imax; i++) 
       for (int j=1; j<=jmax; j++)
           for (int k=1; k<=jmax; k++)
               if (grid_.isFluid(i,j,k))
                  Pmean += p(i,j,k);
   Pmean /= grid_.getNumFluid();
   for (int i=1; i<=imax; i++) 
       for (int j=1; j<=jmax; j++)
           for (int k=1; k<=jmax; k++)
               if (grid_.isFluid(i,j,k))
                  p(i,j,k) -= Pmean;
}

void FluidSimulator3D::computeFGH()
{
   for (int i=1; i<imax; i++) 
       for (int j=1; j<=jmax; j++)
           for (int k=1; k<=jmax; k++)
               f(i,j,k) = grid_.u(i,j,k,CENTER) + dt *  ( Re_1*(d2U_dx2(i,j,k)+d2U_dy2(i,j,k)+d2U_dz2(i,j,k))
                                                         -dU2_dx(i,j,k) - dUV_dy(i,j,k)- dUW_dz(i,j,k) + GX);

   for (int i=1; i<=imax; i++) 
       for (int j=1; j<jmax; j++)
           for (int k=1; k<=jmax; k++)
               g(i,j,k) = grid_.v(i,j,k,CENTER) + dt *  ( Re_1*( d2V_dx2(i,j,k)+d2V_dy2(i,j,k)+d2V_dz2(i,j,k))
                                                         -dV2_dy(i,j,k) - dVU_dx(i,j,k)- dVW_dz(i,j,k) + GY);

   for (int i=1; i<=imax; i++) 
       for (int j=1; j<=jmax; j++)
           for (int k=1; k<jmax; k++)
               h(i,j,k) = grid_.w(i,j,k,CENTER) + dt *  ( Re_1*( d2W_dx2(i,j,k)+d2W_dy2(i,j,k)+d2W_dz2(i,j,k))
                                                         -dW2_dz(i,j,k) - dWU_dx(i,j,k)- dWV_dy(i,j,k) + GZ);

   for (int j=0; j<=jmax; j++)
       for (int k=0; k<=kmax; k++) {
           f(0,j,k) = u(0,j,k);
           f(imax,j,k) = u(imax,j,k); }           
   for (int i=0; i<=imax; i++)
       for (int k=0; k<=kmax; k++) {
           g(i,0,k) = v(i,0,k);
           g(i,jmax,k) = v(i,jmax,k); }        
   for (int i=0; i<=imax; i++)
       for (int j=0; j<=jmax; j++) {
           h(i,j,0) = w(i,j,0);
           h(i,j,kmax) = w(i,j,kmax); }        
}

void FluidSimulator3D::composeRHS()
{
   for (int i=1; i<=imax; i++) 
       for (int j=1; j<=jmax; j++)
           for (int k=1; k<=kmax; k++)
               if (grid_.isFluid(i,j,k))
                  rhs(i,j,k) = (1.0/dt) * ( (grid_.f(i,j,k,CENTER) - grid_.f(i,j,k,WEST))/dx
                                          + (grid_.g(i,j,k,CENTER) - grid_.g(i,j,k,SOUTH))/dy 
                                          + (grid_.h(i,j,k,CENTER) - grid_.h(i,j,k,DOWN))/dz );
}

void FluidSimulator3D::updateVelocities()
{
   for (int i=1; i<imax; i++) 
       for (int j=1; j<=jmax; j++)
           for (int k=1; k<=kmax; k++)
               if (grid_.isFluid(i,j,k) && grid_.isFluid(i+1,j,k))
                  u(i,j,k) = f(i,j,k) - dt/dx *  ( p(i+1,j,k) - p(i,j,k) );
   for (int i=1; i<=imax; i++) 
       for (int j=1; j<jmax; j++)
           for (int k=1; k<=kmax; k++)
               if (grid_.isFluid(i,j,k) && grid_.isFluid(i,j+1,k))
                  v(i,j,k) = g(i,j,k) - dt/dy *  ( p(i,j+1,k) - p(i,j,k) );

   for (int i=1; i<=imax; i++) 
       for (int j=1; j<=jmax; j++)
           for (int k=1; k<kmax; k++)
               if (grid_.isFluid(i,j,k) && grid_.isFluid(i,j,k+1))
                  w(i,j,k) = h(i,j,k) - dt/dz *  ( p(i,j,k+1) - p(i,j,k) );
}

void FluidSimulator3D::simulateTimeStepCount()
{
   VTKWriter3D vtkWriter(grid_, name, true, true);
   while (timeStepNumber<=noOfTimeSteps)
   {
      if (tau>0)
         determineNextDT( dtmax );
      refreshBoundaries();
      if (timeStepNumber%outputInterval==0)
         vtkWriter.write();
      if (timeStepNumber%normFrequency == 0)
         normalizePressure();
      computeFGH();
      composeRHS();
      std::cout<< "\nTime step    = "<<timeStepNumber << "\tTime     = "<< time<<"\n";
      solver_->solve(grid_);   
      updateVelocities();
      timeStepNumber++;
      time += dt;
   }
}

void FluidSimulator3D::simulate( real duration )
{
   VTKWriter3D vtkWriter(grid_, name, true, true);
   while (time<=duration)
   {
      if (tau>0)
         determineNextDT( dtmax );
      refreshBoundaries();
      if (timeStepNumber%outputInterval==0)
         vtkWriter.write();
      if (timeStepNumber%normFrequency == 0)
         normalizePressure();
      computeFGH();
      composeRHS();
      std::cout<< "\nTime step    = "<<timeStepNumber << "\tTime     = "<< time<<"\n";
      solver_->solve(grid_);   
      updateVelocities();
      timeStepNumber++;
      time += dt;
   }
}
