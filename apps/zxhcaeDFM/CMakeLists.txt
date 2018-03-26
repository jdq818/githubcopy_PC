#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDFM zxhcaeDFM.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDFM zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDFM)
