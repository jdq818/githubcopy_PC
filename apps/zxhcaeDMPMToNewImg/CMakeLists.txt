#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMPMToNewImg zxhcaeDMPMToNewImg.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMPMToNewImg zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMPMToNewImg)
