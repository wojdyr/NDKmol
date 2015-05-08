# NOTE: it's only for GLUT (desktop) version of ndkmol
project(NDKmol C CXX)
cmake_minimum_required(VERSION 2.6)

option(USE_ES2 "Use OpenGL ES 2" OFF)
option(USE_NFD "Use nativefiledialog library" ON)
option(DOWNLOAD_ZLIB "Fetch zlib from its website" OFF)
option(DOWNLOAD_FREEGLUT "Fetch freeglut from its website" OFF)
if(UNIX AND NOT APPLE)
 option(USE_X11_ICON "Use icon in X11 Window Manager" ON)
endif()
#if(APPLE)
# option(MACOSX_BUNDLE "Make OSX app bundle" OFF)
#endif()

# TODO: ${THIRDPARTY_INSTALL} is used during 'make' so should be user-writable,
# OTOH libraries are needed at runtime
set(THIRDPARTY_INSTALL "${CMAKE_INSTALL_PREFIX}")

if (POLICY CMP0011) # for FindGTK2
  cmake_policy(SET CMP0011 NEW)
endif()

if ("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC") # avoiding dereferencing MSVC
 add_definitions(-D_USE_MATH_DEFINES)  # needed for M_PI
endif()

include(CheckCXXCompilerFlag)
find_package(OpenGL REQUIRED)

set(LIB_PREFIX ${CMAKE_STATIC_LIBRARY_PREFIX})
set(LIB_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})

set(thirdparty_dir ${CMAKE_SOURCE_DIR}/3rdparty)
if (DOWNLOAD_ZLIB)
  include(ExternalProject)
  ExternalProject_Add(zlib
    URL "http://zlib.net/zlib-1.2.8.tar.gz"
    URL_MD5 44d667c142d7cda120332623eab69f40
    DOWNLOAD_DIR ${thirdparty_dir}
    SOURCE_DIR ${thirdparty_dir}/zlib
    BINARY_DIR ${CMAKE_BINARY_DIR}/zlib
    PREFIX zlib-stamps
    CMAKE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
               -DCMAKE_INSTALL_PREFIX:PATH=${THIRDPARTY_INSTALL}
               -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
               -DSKIP_INSTALL_FILES=ON)
  set(ZLIB_INCLUDE_DIR "${THIRDPARTY_INSTALL}/include")
  if (WIN32)
    set(ZLIB_NAME zlibstatic)
  else()
    set(ZLIB_NAME libz)
  endif()
  set(ZLIB_LIBRARIES "${THIRDPARTY_INSTALL}/lib/${ZLIB_NAME}${LIB_SUFFIX}")
  message(STATUS "Using local zlib: ${ZLIB_LIBRARIES}")
else()
  find_package(ZLIB REQUIRED)
endif()

if (DOWNLOAD_FREEGLUT)
  include(ExternalProject)
  ExternalProject_Add(freeglut
    URL "http://downloads.sourceforge.net/project/freeglut/freeglut/3.0.0/freeglut-3.0.0.tar.gz"
    URL_MD5 90c3ca4dd9d51cf32276bc5344ec9754
    DOWNLOAD_DIR ${thirdparty_dir}
    SOURCE_DIR ${thirdparty_dir}/freeglut
    BINARY_DIR ${CMAKE_BINARY_DIR}/freeglut
    PREFIX freeglut-stamps
    CMAKE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
               -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
               -DCMAKE_INSTALL_PREFIX:PATH=${THIRDPARTY_INSTALL}
               -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
               -DCMAKE_INSTALL_LIBDIR:PATH=${THIRDPARTY_INSTALL}/lib
               -DFREEGLUT_REPLACE_GLUT=OFF
               -DFREEGLUT_BUILD_SHARED_LIBS=OFF
               -DINSTALL_PDB=OFF)
  set(GLUT_INCLUDE_DIR "${THIRDPARTY_INSTALL}/include")
  set(GLUT_LIBRARIES
      "${THIRDPARTY_INSTALL}/lib/${LIB_PREFIX}freeglut${LIB_SUFFIX}")
  if(UNIX AND NOT(ANDROID OR BLACKBERRY))
    find_package(X11 REQUIRED)
    foreach (_lib Xrandr Xxf86vm Xinput)
      if(X11_${_lib}_LIB)
        list(APPEND GLUT_LIBRARIES ${X11_${_lib}_LIB})
      endif()
    endforeach()
    list(APPEND GLUT_LIBRARIES ${X11_LIBRARIES})
  endif()
  message(STATUS "Using local freeglut: ${GLUT_LIBRARIES}")
else()
  find_package(GLUT REQUIRED)
  # FindGLUT.cmake adds libXmu and libXi even if it's not needed
  if (NOT GLUT_Xmu_LIBRARY)
    set(GLUT_Xmu_LIBRARY "")
  endif()
  if (NOT GLUT_Xi_LIBRARY)
    set(GLUT_Xi_LIBRARY "")
  endif()
  set(GLUT_LIBRARIES ${GLUT_glut_LIBRARY} ${GLUT_Xmu_LIBRARY} ${GLUT_Xi_LIBRARY}
                     ${GLUT_cocoa_LIBRARY})
endif()

if (DEFINED USE_X11_ICON AND USE_X11_ICON)
  find_package(X11 REQUIRED)
  if (NOT X11_Xpm_LIB)
    message(FATAL_ERROR "libXpm not found. On Linux you need a package named "
                        "libxpm-dev or libXpm-devel or similarly.")
  endif()
  add_definitions(-DUSE_X11_ICON)
  set(EXTRA_LIBS ${X11_Xpm_LIB})
endif()

if (USE_NFD)
  add_definitions(-DUSE_NFD)
  find_path(NFD_INCLUDE nfd.h)
  if (NOT NFD_INCLUDE)
    include(ExternalProject)
    ExternalProject_Add(nativefiledialog
      GIT_REPOSITORY https://github.com/mlabbe/nativefiledialog
      SOURCE_DIR ${CMAKE_SOURCE_DIR}/nfd
      PREFIX nfd-stamps
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND "")
    set(NFD_SRC ${CMAKE_SOURCE_DIR}/nfd/src)
    set(NFD_INCLUDE ${NFD_SRC}/include)
    include_directories(${NFD_INCLUDE})
    if (WIN32)
      set(NFD_PLAT_SRC ${NFD_SRC}/nfd_win.cpp)
    elseif (APPLE)
      set(NFD_PLAT_SRC ${NFD_SRC}/nfd_cocoa.m)
    else()
      find_package(GTK2 2.6 REQUIRED gtk)
      set(NFD_PLAT_SRC ${NFD_SRC}/nfd_gtk.c)
    endif()
    set_source_files_properties(${NFD_SRC}/nfd_common.c ${NFD_PLAT_SRC}
                                PROPERTIES GENERATED TRUE)
    add_library(nfd STATIC ${NFD_SRC}/nfd_common.c ${NFD_PLAT_SRC})
    # neither set_property nor set_target_properties seem to work with
    # cmake 2.8.7, so we use include_directories instead
    #set_property(TARGET nfd PROPERTY INCLUDE_DIRECTORIES
    #             ${NFD_INCLUDE} ${GTK2_INCLUDE_DIRS})
    include_directories(${GTK2_INCLUDE_DIRS})
    set(EXTRA_LIBS ${EXTRA_LIBS} nfd ${GTK2_LIBRARIES})
  endif()
endif()

include_directories(${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIR})
include_directories(${ZLIB_INCLUDE_DIRS} ${X11_INCLUDE_DIR})
include_directories("${PROJECT_SOURCE_DIR}")

if (WIN32 AND NOT EXISTS "${CMAKE_SOURCE_DIR}/GL/glext.h")
  file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/GL")
  message(STATUS "Downloading glext.h...")
  file(DOWNLOAD https://www.opengl.org/registry/api/GL/glext.h
       "${CMAKE_SOURCE_DIR}/GL/glext.h" SHOW_PROGRESS)
endif()

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

add_definitions("-DINITIAL_PDB=\"${CMAKE_INSTALL_PREFIX}/share/ndkmol/initial.pdb\"")

if(WIN32)
  set(RCFILE GLUT/ndkmol.rc)
endif()

add_executable(ndkmol WIN32
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
GLUT/NDKmol.cpp
${RCFILE}
res/raw/initial.pdb
)

if (USE_NFD)
  add_dependencies(ndkmol nativefiledialog)
endif()
if (DOWNLOAD_ZLIB)
  add_dependencies(ndkmol zlib)
endif()
if (DOWNLOAD_FREEGLUT)
  add_dependencies(ndkmol freeglut)
endif()

# We could make an app bundle, but would it be useful?
#if (APPLE AND MACOSX_BUNDLE)
#  set_target_properties(ndkmol PROPERTIES
#                        MACOSX_BUNDLE TRUE
#                        MACOSX_BUNDLE_BUNDLE_NAME "NDKmol")
#  set(MACOSX_BUNDLE_ICON_FILE ndkmol.icns)
#  set_source_files_properties(res/raw/initial.pdb PROPERTIES
#                              MACOSX_PACKAGE_LOCATION Resources)
#endif()

target_link_libraries(ndkmol ${EXTRA_LIBS} ${GLUT_LIBRARIES}
                             ${OPENGL_LIBRARIES} ${ZLIB_LIBRARIES})

# check if the program starts at least
add_custom_target(check COMMAND ndkmol "--version")

install(TARGETS ndkmol DESTINATION bin)
install(FILES res/raw/initial.pdb DESTINATION share/ndkmol/)
