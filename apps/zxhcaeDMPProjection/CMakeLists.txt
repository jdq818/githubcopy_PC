#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMPProjection zxhcaeDMPProjection.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMPProjection zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMPProjection)
