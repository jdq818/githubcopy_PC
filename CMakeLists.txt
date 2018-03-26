cmake_minimum_required(VERSION 2.6)
#
# Set the project name.
#

PROJECT (zxhprojCAE)



SET (zxhprojCAE_CAN_BUILD 1)

SET (EXECUTABLE_OUTPUT_PATH  ${PROJECT_BINARY_DIR}/bin CACHE INTERNAL
  "Where to put the executables")

SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib CACHE INTERNAL
  "Where to put the libraries")

 

#-----------------------------------------------------------------------------
# Setup the path to load our macros
#-----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH "${UCLTK_BINARY_DIR}/CMakeModules")
   

IF (zxhprojCAE_CAN_BUILD)

  #
  # Here is where you can list the sub-directories holding your local 
  # classes. Sorting classes by 'package' type like VTK does (Common, 
  # Rendering, Filtering, Imaging, IO, etc.) is a good thing and prevents 
  # numerous dependencies problems.
  #

  SUBDIRS ( 
  class 
  ExternLibs
  apps
  )
 
  #
  # You can put your include path(s) here
  #

  
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/class)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/ExternLibs) 
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})
  
  
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/ExternLibs/nifti/niftilib)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/ExternLibs/nifti/znzlib)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/ExternLibs/Common)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/ExternLibs/Images)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/ExternLibs/ImageModelings)


INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/ExternLibs/nifti/zlib) 
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/ExternLibs/nifti/zlib) 

FIND_PACKAGE(VTK)
INCLUDE(${VTK_USE_FILE})

#INCLUDE_DIRECTORIES(${VTK_DIR}/include/vtk-5.10)

#LINK_DIRECTORIES(${VTK_DIR}/lib/vtk-5.10-dbg)

add_definitions(-DHAVE_NIFTI)
add_definitions(-DHAVE_ZLIB)
 
  
ENDIF (zxhprojCAE_CAN_BUILD)
