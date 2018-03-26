#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMP_1.1_MLLMLR zxhcaeDMP_1.1_MLLMLR.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMP_1.1_MLLMLR zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMP_1.1_MLLMLR)
