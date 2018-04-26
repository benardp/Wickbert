project( Wickbert )


cmake_minimum_required(VERSION 2.8.11)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#-----------------------------------------------------------------------------
# Test for some required system information.
#INCLUDE (${CMAKE_ROOT}/Modules/CMakeBackwardCompatibilityC.cmake)
#INCLUDE (${CMAKE_ROOT}/Modules/CMakeBackwardCompatibilityCXX.cmake)

#-----------------------------------------------------------------------------
# Output directories.
SET (LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin CACHE INTERNAL "Single output directory for building all libraries.")
SET (EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin CACHE INTERNAL "Single output directory for building all executables.")
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)
SET(GPS_LIBRARY_PATH "${LIBRARY_OUTPUT_PATH}")
SET(GPS_EXECUTABLE_PATH "${EXECUTABLE_OUTPUT_PATH}")
SET(CXX_TEST_PATH ${EXECUTABLE_OUTPUT_PATH})

#Create necessary directories and files 
FILE(MAKE_DIRECTORY  ${PROJECT_SOURCE_DIR}/bin/release )
FILE(MAKE_DIRECTORY  ${PROJECT_SOURCE_DIR}/bin/debug )
FILE(MAKE_DIRECTORY  ${PROJECT_SOURCE_DIR}/bin/release/bitmaps )
FILE(MAKE_DIRECTORY  ${PROJECT_SOURCE_DIR}/bin/debug/bitmaps )

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/wxModeler/defaultparticlesystem.par  ${PROJECT_SOURCE_DIR}/bin/release/defaultparticlesystem.par  COPYONLY)
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/wxModeler/defaultparticlesystem.par  ${PROJECT_SOURCE_DIR}/bin/debug/defaultparticlesystem.par  COPYONLY)

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/wxModeler/bitmaps/Particles.bmp  ${PROJECT_SOURCE_DIR}/bin/release/bitmaps/Particles.bmp  COPYONLY)
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/wxModeler/bitmaps/Particles.bmp  ${PROJECT_SOURCE_DIR}/bin/debug/bitmaps/Particles.bmp  COPYONLY)

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/wxModeler/bitmaps/Particles.xpm  ${PROJECT_SOURCE_DIR}/bin/release/bitmaps/Particles.xpm  COPYONLY)
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/wxModeler/bitmaps/Particles.xpm  ${PROJECT_SOURCE_DIR}/bin/debug/bitmaps/Particles.xpm  COPYONLY)

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/wxModeler/bitmaps/ParticleStuff.bmp  ${PROJECT_SOURCE_DIR}/bin/release/bitmaps/ParticleStuff.bmp  COPYONLY)
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/wxModeler/bitmaps/ParticleStuff.bmp  ${PROJECT_SOURCE_DIR}/bin/debug/bitmaps/ParticleStuff.bmp  COPYONLY)

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/wxModeler/bitmaps/ParticleStuffHeader.bmp  ${PROJECT_SOURCE_DIR}/bin/release/bitmaps/ParticleStuffHeader.bmp  COPYONLY)
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/wxModeler/bitmaps/ParticleStuffHeader.bmp  ${PROJECT_SOURCE_DIR}/bin/debug/bitmaps/ParticleStuffHeader.bmp  COPYONLY)

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/wxModeler/bitmaps/ParticleSystem.bmp  ${PROJECT_SOURCE_DIR}/bin/release/bitmaps/ParticleSystem.bmp  COPYONLY)
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/wxModeler/bitmaps/ParticleSystem.bmp  ${PROJECT_SOURCE_DIR}/bin/debug/bitmaps/ParticleSystem.bmp  COPYONLY)

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/wxModeler/bitmaps/ParticleSystem.xpm  ${PROJECT_SOURCE_DIR}/bin/release/bitmaps/ParticleSystem.xpm  COPYONLY)
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/wxModeler/bitmaps/ParticleSystem.xpm  ${PROJECT_SOURCE_DIR}/bin/debug/bitmaps/ParticleSystem.xpm  COPYONLY)

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/wxModeler/bitmaps/ParticleSystems.bmp  ${PROJECT_SOURCE_DIR}/bin/release/bitmaps/ParticleSystems.bmp  COPYONLY)
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/wxModeler/bitmaps/ParticleSystems.bmp  ${PROJECT_SOURCE_DIR}/bin/debug/bitmaps/ParticleSystems.bmp  COPYONLY)

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/wxModeler/bitmaps/ParticleSystems.xpm  ${PROJECT_SOURCE_DIR}/bin/release/bitmaps/ParticleSystems.xpm  COPYONLY)
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/wxModeler/bitmaps/ParticleSystems.xpm  ${PROJECT_SOURCE_DIR}/bin/debug/bitmaps/ParticleSystems.xpm  COPYONLY)

#-----------------------------------------------------------------------------
# GPS build configuration options.

OPTION(WB_USE_VTK "Some Implcits use VTK. Set it ON if you have VTK installed" OFF)
OPTION(WB_USE_ITK "Some Implcits use ITK. Set it ON if you have ITK installed" OFF)
OPTION(WB_USE_CGAL "Some Shaders use CGAL. Set it ON if you have CGAL installed" OFF)
OPTION(WB_USE_CG "Some Shaders use NVIDIA CG (WIN32 only). Set it ON if you have CG installed" OFF)
OPTION(WB_USE_SFL "MLS can use SFL. Set it ON if you have SFL installed" OFF)

IF(WIN32)
   IF(BUILD_SHARED_LIBS)
     SET(GPS_COMMON_BUILD_TYPE "SHARED")
   ELSE(BUILD_SHARED_LIBS)
     SET(GPS_COMMON_BUILD_TYPE "STATIC")
   ENDIF(BUILD_SHARED_LIBS)
   SET(BUILD_SHARED_LIBS OFF)
ENDIF(WIN32)

#Find packages 
FIND_PACKAGE(OpenGL)

IF(WB_USE_ITK)
    FIND_PACKAGE(ITK)
    IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
  ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
  "ITK not found. Please set ITK_DIR.")
    ENDIF(ITK_FOUND)
ENDIF(WB_USE_ITK)

IF(WB_USE_VTK)
    FIND_PACKAGE(VTK)
    IF(VTK_FOUND)
      INCLUDE(${VTK_USE_FILE})
  ELSE(VTK_FOUND)
      MESSAGE(FATAL_ERROR
          "VTK not found. Please set VTK_DIR.")
    ENDIF(VTK_FOUND)
ENDIF(WB_USE_VTK)

# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()
endif()

FIND_PACKAGE(GSL)

IF(WB_USE_SFL)
    SET(SFL_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/sfl/sfl" CACHE PATH "Path to SFL include files" )
    SET(SFL_LIB_PATH "${PROJECT_SOURCE_DIR}/lib-release" CACHE PATH "Path to SFL lib" )
    SET(SFL_LIB_DEBUG_PATH "${PROJECT_SOURCE_DIR}/lib-debug" CACHE PATH "Path to Debug SFL lib" )
    SET(PI_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/sfl/common/pi" CACHE PATH "Path to PI include files" )
ENDIF(WB_USE_SFL)

IF( WIN32)
    SET(GPS_RC_FILE wxModeler.rc)
ELSE( WIN32)
    SET(GPS_RC_FILE "")
ENDIF( WIN32)

#CMake Options
SET(CMAKE_VERBOSE_MAKEFILE TRUE)

find_package(wxWidgets 3.1 REQUIRED gl core base)
include(${wxWidgets_USE_FILE})

find_package(OpenMesh REQUIRED)

IF(WB_USE_CG)
SET(CG_LIB_PATH CACHE PATH "Path to NVIDIA CG library" )
SET(CG_INCLUDE_PATH CACHE PATH "Path to NVIDIA CG includes" )
SET(GLH_INCLUDE_PATH CACHE PATH "Path to nvsdk GLH includes" )
INCLUDE_DIRECTORIES(
  ${CG_INCLUDE_PATH} ${GLH_INCLUDE_PATH}
)
ENDIF(WB_USE_CG)

IF(WB_USE_CGAL)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )
if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  
endif()

# include helper file
include( ${CGAL_USE_FILE} )

# Creating entries for all C++ files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )
ENDIF(WB_USE_CGAL)

set(GPS_SOURCE_DIR  ${PROJECT_SOURCE_DIR})

#-----------------------------------------------------------------------------
# Dispatch the build into the proper subdirectories.
SUBDIRS(ism libgm Particles Surface wxModeler tests)
INCLUDE_DIRECTORIES(
  ism 
  libgm 
  Particles 
  Surface
  ${PROJECT_SOURCE_DIR}
  ${OPENGL_INCLUDE_DIR}
  ${WX_H_INCLUDE_PATH} 
  ${WX_CONTRIB_INCLUDE_PATH}
  ${WX_INCLUDE_PATH} 
  ${Boost_INCLUDE_DIRS}
  ${OPENMESH_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIRS}
  ${SFL_INCLUDE_PATH}
  ${PI_INCLUDE_PATH}  
  ${SFL_INCLUDE_PATH}/pi
)

#-----------------------------------------------------------------------------
# GPS requires special compiler flags on some platforms.
IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(GPS_REQUIRED_CXX_FLAGS "${GPS_REQUIRED_CXX_FLAGS} -ftemplate-depth-50")
  IF(APPLE)
    SET(GPS_REQUIRED_C_FLAGS "${GPS_REQUIRED_C_FLAGS} -no-cpp-precomp")
    SET(GPS_REQUIRED_CXX_FLAGS "${GPS_REQUIRED_CXX_FLAGS} -no-cpp-precomp")
    SET(GPS_REQUIRED_CXX_FLAGS "${GPS_REQUIRED_CXX_FLAGS} -Wno-long-double")
    SET(GPS_REQUIRED_CXX_FLAGS "${GPS_REQUIRED_CXX_FLAGS} -Wno-deprecated")  
  ENDIF(APPLE)
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

IF(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
ENDIF(UNIX)

#wxwidgets requires  D_USE_MATH_DEFINES for WIN32 
IF(WIN32)
  ADD_DEFINITIONS(-D_USE_MATH_DEFINES -DNOMINMAX -D_WIN32)
ENDIF(WIN32)

#VS2005 is picky
ADD_DEFINITIONS(-D_HAS_ITERATOR_DEBUGGING=0
        -D_SECURE_SCL=0 
        -D_CRT_SECURE_NO_DEPRECATE
        )
IF(WB_USE_ITK)
  ADD_DEFINITIONS(-DWB_USE_ITK)
ENDIF(WB_USE_ITK)

IF(WB_USE_VTK)
  ADD_DEFINITIONS(-DWB_USE_VTK)
ENDIF(WB_USE_VTK)

IF(WB_USE_CGAL)
  ADD_DEFINITIONS(-DWB_USE_CGAL)
ENDIF(WB_USE_CGAL)

IF(WB_USE_CG)
  ADD_DEFINITIONS(-DWB_USE_CG)
ENDIF(WB_USE_CG)

IF(WB_USE_SFL)
  ADD_DEFINITIONS(-DWB_USE_SFL)
ENDIF(WB_USE_SFL)

LINK_DIRECTORIES(${WX_LIB_PATH} ${CG_LIB_PATH} ${OPENMESH_LIBRARY_DIR} ${CGAL_LIB_PATH} ${SFL_LIB_PATH} ${GMP_LIB_PATH})

IF(WIN32)
LINK_DIRECTORIES(${OPENMESH_LIB_DEBUG_PATH} ${GSL_LIB_DEBUG_PATH} 
${SFL_LIB_DEBUG_PATH})
ENDIF(WIN32)

IF(UNIX AND NOT APPLE)
  ADD_DEFINITIONS(-DGTK_NO_CHECK_CASTS -D__WXGTK__ -D_FILE_OFFSET_BITS=64 -D_LARGE_FILES 
-D_LARGEFILE_SOURCE=1 -DNO_GCC_PRAGMA -DINCLUDE_TEMPLATES -DARCH_LINUX)
ENDIF(UNIX AND NOT APPLE)


