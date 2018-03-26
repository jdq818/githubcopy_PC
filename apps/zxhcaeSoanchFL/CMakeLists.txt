#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeSoanchFL zxhcaeSoanchFL.cpp)
TARGET_LINK_LIBRARIES(zxhcaeSoanchFL zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeSoanchFL)
