#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMPResoImageU1 zxhcaeDMPResoImageU1.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMPResoImageU1 zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMPResoImageU1)
