#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMPCMBORIIMG zxhcaeDMPCMBORIIMG.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMPCMBORIIMG zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMPCMBORIIMG)
