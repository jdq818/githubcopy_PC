#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDFM_MUS zxhcaeDFM_MUS.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDFM_MUS zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDFM_MUS)
