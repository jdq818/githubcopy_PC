#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMP_1.0_WMMLL zxhcaeDMP_1.0_WMMLL.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMP_1.0_WMMLL zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMP_1.0_WMMLL)
