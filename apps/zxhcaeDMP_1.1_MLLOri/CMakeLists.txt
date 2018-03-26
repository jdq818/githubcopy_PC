#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMP_1.1_MLLOri zxhcaeDMP_1.1_MLLOri.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMP_1.1_MLLOri zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMP_1.1_MLLOri)
