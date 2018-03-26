#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMPSmooth zxhcaeDMPSmooth.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMPSmooth zxhImages zxhCommon znz zlib niftiio vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMPSmooth)
