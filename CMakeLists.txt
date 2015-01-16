# NOTE: it's only for GLUT (desktop) version of ndkmol
cmake_minimum_required(VERSION 2.8)

include(CheckCXXCompilerFlag)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(ZLIB REQUIRED)
include_directories(${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS})
include_directories(${ZLIB_INCLUDE_DIRS})
include_directories("${PROJECT_SOURCE_DIR}")

option(USE_ES2 "Use OpenGL ES 2" ON)

# show warnings by default
check_cxx_compiler_flag(-Wextra HAS_WEXTRA_FLAG)
if (HAS_WEXTRA_FLAG)
  set(EXTRA_CXX_FLAGS "-Wall -Wextra" CACHE STRING "")
endif()
set(EXTRA_CXX_FLAGS ${EXTRA_CXX_FLAGS} CACHE STRING "Flags for compiler" FORCE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_CXX_FLAGS}")

# set_property(TARGET ndkmol PROPERTY CXX_STANDARD 11) # only in cmake 3.1+
check_cxx_compiler_flag(-std=c++11 HAS_STDCXX11_FLAG)
if (HAS_STDCXX11_FLAG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

if (NOT USE_ES2)
  add_definitions(-DOPENGL_ES1)
endif()

add_executable(ndkmol
NDKmol/Atom.cpp              NDKmol/PDBReader.cpp
NDKmol/CCP4Reader.cpp        NDKmol/Protein.cpp
NDKmol/ChemDatabase.cpp      NDKmol/Quaternion.cpp
NDKmol/Color.cpp             NDKmol/Renderable.cpp
NDKmol/CylinderGeometry.cpp  NDKmol/RibbonStrip.cpp
NDKmol/Debug.cpp             NDKmol/SDFReader.cpp
NDKmol/Geometry.cpp          NDKmol/SmoothCurve.cpp
NDKmol/GLES.cpp              NDKmol/SmoothTube.cpp
NDKmol/Line.cpp              NDKmol/SphereGeometry.cpp
NDKmol/MarchingSquares.cpp   NDKmol/VBOCylinder.cpp
NDKmol/MatRenderable.cpp     NDKmol/VBOSphere.cpp
NDKmol/Matrix.cpp            NDKmol/Vector3.cpp
NDKmol/NdkView.cpp           #NDKmol/VolumeRenderer.cpp
GLUT/NDKmol.cpp)

target_link_libraries(ndkmol ${GLUT_LIBRARY} ${OPENGL_LIBRARIES}
                             ${ZLIB_LIBRARIES})

install(TARGETS ndkmol DESTINATION bin)
