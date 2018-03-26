#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMP_1.0_WMMLLT zxhcaeDMP_1.0_WMMLLT.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMP_1.0_WMMLLT zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMP_1.0_WMMLLT)
