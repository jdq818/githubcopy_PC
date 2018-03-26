#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMPCalcDistFREF zxhcaeDMPCalcDistFREF.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMPCalcDistFREF zxhImages zxhCommon znz zlib niftiio vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMPCalcDistFREF)
