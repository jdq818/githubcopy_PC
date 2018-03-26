#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDFMPO zxhcaeDFMPO.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDFMPO zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDFMPO)
