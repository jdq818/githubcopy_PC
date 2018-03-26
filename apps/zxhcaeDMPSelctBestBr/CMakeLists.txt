#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMPSelctBestBr zxhcaeDMPSelctBestBr.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMPSelctBestBr zxhImages zxhCommon znz zlib niftiio vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMPSelctBestBr)
