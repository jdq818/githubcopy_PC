#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMP_1.1_MLLCOWM zxhcaeDMP_1.1_MLLCOWM.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMP_1.1_MLLCOWM zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMP_1.1_MLLCOWM)
