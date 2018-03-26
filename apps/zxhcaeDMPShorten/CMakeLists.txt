#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMPShorten zxhcaeDMPShorten.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMPShorten zxhImages zxhCommon znz zlib niftiio vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMPShorten)
