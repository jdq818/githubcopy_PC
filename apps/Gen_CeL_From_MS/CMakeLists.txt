#
# Add the executable
#

ADD_EXECUTABLE(Gen_CeL_From_MS Gen_CeL_From_MS.cpp)
TARGET_LINK_LIBRARIES(Gen_CeL_From_MS zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin Gen_CeL_From_MS)
