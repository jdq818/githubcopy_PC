#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMPSelctBestBr4T zxhcaeDMPSelctBestBr4T.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMPSelctBestBr4T zxhImages zxhCommon znz zlib niftiio vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMPSelctBestBr4T)
