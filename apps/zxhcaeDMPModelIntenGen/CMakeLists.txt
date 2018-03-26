#
# Add the executable
#

ADD_EXECUTABLE(zxhcaeDMPModelIntenGen zxhcaeDMPModelIntenGen.cpp)
TARGET_LINK_LIBRARIES(zxhcaeDMPModelIntenGen zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin zxhcaeDMPModelIntenGen)
