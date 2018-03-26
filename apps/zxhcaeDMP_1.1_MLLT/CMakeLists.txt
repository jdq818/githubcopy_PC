#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMP_1.1_MLLT zxhcaeDMP_1.1_MLLT.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMP_1.1_MLLT zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMP_1.1_MLLT)
