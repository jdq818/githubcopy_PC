#
# Add the executable
#

ADD_EXECUTABLE(zxhcaetransImagetool zxhcaetransImagetool.cpp)
TARGET_LINK_LIBRARIES(zxhcaetransImagetool zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaetransImagetool)
