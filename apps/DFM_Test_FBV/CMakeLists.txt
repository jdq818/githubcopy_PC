#
# Add the executable
#

ADD_EXECUTABLE(DFM_Test_FBV DFM_Test_FBV.cpp)
TARGET_LINK_LIBRARIES(DFM_Test_FBV zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin DFM_Test_FBV)
