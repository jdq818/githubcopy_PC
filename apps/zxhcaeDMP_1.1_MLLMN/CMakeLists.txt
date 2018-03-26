#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMP_1.1_MLLMN zxhcaeDMP_1.1_MLLMN.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMP_1.1_MLLMN zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMP_1.1_MLLMN)
