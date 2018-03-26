#
# Add the executable
#

ADD_EXECUTABLE(miiMinPathModelApp miiMinPathModelApp.cpp)
TARGET_LINK_LIBRARIES(miiMinPathModelApp zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin miiMinPathModelApp)
