#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDFM_KT zxhcaeDFM_KT.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDFM_KT zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDFM_KT)
