#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMPMultiModel zxhcaeDMPMultiModel.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMPMultiModel zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMPMultiModel)
