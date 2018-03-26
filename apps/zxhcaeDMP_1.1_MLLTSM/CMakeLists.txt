#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMP_1.1_MLLTSM zxhcaeDMP_1.1_MLLTSM.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMP_1.1_MLLTSM zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMP_1.1_MLLTSM)
