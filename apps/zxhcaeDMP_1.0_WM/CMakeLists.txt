#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMP_1.0_WM zxhcaeDMP_1.0_WM.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMP_1.0_WM zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMP_1.0_WM)
