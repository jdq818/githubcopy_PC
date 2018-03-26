#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMPBesselCurve zxhcaeDMPBesselCurve.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMPBesselCurve zxhImages zxhCommon znz zlib niftiio vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMPBesselCurve)
