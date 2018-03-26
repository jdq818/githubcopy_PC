#
# Add the executable
#

ADD_EXECUTABLE(Curfusion_jdq2018_v1 Curfusion_jdq2018_v1.cpp)
TARGET_LINK_LIBRARIES(Curfusion_jdq2018_v1 zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin Curfusion_jdq2018_v1)
