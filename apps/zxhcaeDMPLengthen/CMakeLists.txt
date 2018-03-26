#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMPLengthen zxhcaeDMPLengthen.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMPLengthen zxhImages zxhCommon znz zlib niftiio vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMPLengthen)
