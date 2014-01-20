#ifndef VTKFILEWRITER_HH
#define VTKFILEWRITER_HH


#include "StaggeredGrid3D.hh"



//*******************************************************************************************************************
/*! Writes StaggeredGrid as vtk file

  - vtk files can for example be opened with Paraview ( http://www.paraview.org/ )
  - writes out pressure and/or velocity

  - Usage:
   \code
       VTKWriter vtkWriter ( myGrid, "lidDrivenCavity", true, true );
       // for each timestep:
      vtkWriter.write();
    \endcode
    This creates on file per timestep: "lidDrivenCavity_0001.vtk", ""lidDrivenCavity_0002.vtk" ...

*/
//*******************************************************************************************************************
class VTKWriter3D
{

public:

   VTKWriter3D(  const StaggeredGrid3D & grid, const std::string & basename,
               bool writePressure = true, bool writeVelocity = true );

   void write();

private:
   const StaggeredGrid3D & grid_;
   std::string baseName_;

   bool writeVelocity_;
   bool writePressure_;

   int counter_;
   std::string header_;

};



#endif




