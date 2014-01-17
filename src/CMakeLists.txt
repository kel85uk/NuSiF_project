
include_directories( . )

add_executable( nusif FileReader.cc Debug.cc Array.hh lodepng.cc GrayScaleImage.cc StaggeredGrid.cc SORSolver.cc RedBlackSORSolver.cc CGSolver.cc MultigridSolver.cc FluidSimulator.cc VTKWriter.cc nusif.cc)
add_executable( nusif3D FileReader.cc Debug.cc Array.hh StaggeredGrid3D.cc SORSolver3D.cc FluidSimulator3D.cc VTKWriter3D.cc nusif3D.cc)
