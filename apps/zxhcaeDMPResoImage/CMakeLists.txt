#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMPResoImage zxhcaeDMPResoImage.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMPResoImage zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMPResoImage)
