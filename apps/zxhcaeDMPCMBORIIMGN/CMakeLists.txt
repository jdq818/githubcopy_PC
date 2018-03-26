#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMPCMBORIIMGN zxhcaeDMPCMBORIIMGN.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMPCMBORIIMGN zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMPCMBORIIMGN)
