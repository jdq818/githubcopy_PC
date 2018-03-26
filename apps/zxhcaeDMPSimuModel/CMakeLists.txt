#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMPSimuModel zxhcaeDMPSimuModel.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMPSimuModel zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMPSimuModel)
