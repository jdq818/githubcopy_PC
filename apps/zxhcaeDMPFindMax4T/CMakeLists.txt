#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMPFindMax4T zxhcaeDMPFindMax4T.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMPFindMax4T zxhImages zxhCommon znz zlib niftiio vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMPFindMax4T)
