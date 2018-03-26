#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMPMultiResolution zxhcaeDMPMultiResolution.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMPMultiResolution zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMPMultiResolution)
