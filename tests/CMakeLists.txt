include_directories( ../src )
add_executable( Array_test      ../src/Array.hh      ../src/Debug.cc ArrayTest.cc      )

add_executable( Filereader_test ../src/FileReader.cc ../src/Debug.cc FileReaderTest.cc )

add_executable( Derivative_test         ../src/FileReader.cc ../src/Debug.cc ../src/Array.hh  ../src/lodepng.cc ../src/GrayScaleImage.cc ../src/StaggeredGrid.cc ../src/SORSolver.cc ../src/RedBlackSORSolver.cc ../src/CGSolver.cc ../src/MultigridSolver.cc ../src/VTKWriter.cc ../src/FluidSimulator.cc DerivativeTest.cc)

add_executable( SOR_test         ../src/FileReader.cc ../src/Debug.cc ../src/Array.hh ../src/lodepng.cc ../src/GrayScaleImage.cc ../src/StaggeredGrid.cc ../src/SORSolver.cc SORTest.cc)

add_executable( RedBlackSOR_test         ../src/FileReader.cc ../src/Debug.cc ../src/Array.hh ../src/lodepng.cc ../src/GrayScaleImage.cc ../src/StaggeredGrid.cc ../src/RedBlackSORSolver.cc RedBlackSORTest.cc)

add_executable( CG_test         ../src/FileReader.cc ../src/Debug.cc ../src/Array.hh ../src/lodepng.cc ../src/GrayScaleImage.cc ../src/StaggeredGrid.cc ../src/CGSolver.cc CGTest.cc)

add_executable( PCG_test         ../src/FileReader.cc ../src/Debug.cc ../src/Array.hh ../src/MatrixCOO.cc ../src/double_compares.cpp ../src/lodepng.cc ../src/GrayScaleImage.cc ../src/StaggeredGrid.cc ../src/PCGSolver.cc PCGTest.cc)

add_executable( Multigrid_test         ../src/FileReader.cc ../src/Debug.cc ../src/Array.hh ../src/lodepng.cc ../src/GrayScaleImage.cc ../src/StaggeredGrid.cc ../src/MultigridSolver.cc MultigridTest.cc)

add_executable( Derivative_test3D         ../src/FileReader.cc ../src/Debug.cc ../src/Array.hh ../src/StaggeredGrid3D.cc ../src/SORSolver3D.cc ../src/RedBlackSORSolver3D.cc ../src/CGSolver3D.cc ../src/VTKWriter3D.cc ../src/FluidSimulator3D.cc DerivativeTest3D.cc)

add_executable( SOR_test3D         ../src/FileReader.cc ../src/Debug.cc ../src/Array.hh ../src/StaggeredGrid3D.cc ../src/SORSolver3D.cc SORTest3D.cc)

add_executable( RedBlackSOR_test3D         ../src/FileReader.cc ../src/Debug.cc ../src/Array.hh ../src/StaggeredGrid3D.cc ../src/RedBlackSORSolver3D.cc RedBlackSORTest3D.cc)

add_executable( CG_test3D         ../src/FileReader.cc ../src/Debug.cc ../src/Array.hh ../src/StaggeredGrid3D.cc ../src/CGSolver3D.cc CGTest3D.cc)


