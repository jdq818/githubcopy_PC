#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMP_1.1_MLL zxhcaeDMP_1.1_MLL.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMP_1.1_MLL zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMP_1.1_MLL)
