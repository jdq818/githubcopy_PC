#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMPResoImageN zxhcaeDMPResoImageN.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMPResoImageN zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMPResoImageN)
