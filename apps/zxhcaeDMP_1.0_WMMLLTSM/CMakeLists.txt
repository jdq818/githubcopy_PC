#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMP_1.0_WMMLLTSM zxhcaeDMP_1.0_WMMLLTSM.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMP_1.0_WMMLLTSM zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMP_1.0_WMMLLTSM)
