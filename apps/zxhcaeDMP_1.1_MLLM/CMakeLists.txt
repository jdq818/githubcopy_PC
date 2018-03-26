#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMP_1.1_MLLM zxhcaeDMP_1.1_MLLM.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMP_1.1_MLLM zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMP_1.1_MLLM)
