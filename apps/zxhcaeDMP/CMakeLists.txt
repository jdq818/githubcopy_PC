#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMP zxhcaeDMP.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMP zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMP)
