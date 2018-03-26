#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeSoanch zxhcaeSoanch.cpp)
TARGET_LINK_LIBRARIES(zxhcaeSoanch zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeSoanch)
