#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMP_1.1 zxhcaeDMP_1.1.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMP_1.1 zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMP_1.1)
