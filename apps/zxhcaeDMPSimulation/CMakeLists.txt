#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMPSimulation zxhcaeDMPSimulation.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMPSimulation zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMPSimulation)
