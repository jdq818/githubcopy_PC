#
# Add the executable
#

ADD_EXECUTABLE(jdqLineop_resam jdqLineop_resam.cpp)
TARGET_LINK_LIBRARIES(jdqLineop_resam zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin jdqLineop_resam)
