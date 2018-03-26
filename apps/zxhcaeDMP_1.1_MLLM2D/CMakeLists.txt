#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMP_1.1_MLLM2D zxhcaeDMP_1.1_MLLM2D.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMP_1.1_MLLM2D zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMP_1.1_MLLM2D)
