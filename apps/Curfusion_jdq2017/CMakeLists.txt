#
# Add the executable
#

ADD_EXECUTABLE(Curfusion_jdq2017 Curfusion_jdq2017.cpp)
TARGET_LINK_LIBRARIES(Curfusion_jdq2017 zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin Curfusion_jdq2017)
