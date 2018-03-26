#
# Add the executable
#

ADD_EXECUTABLE(JDQ_Curve_Fusion JDQ_Curve_Fusion.cpp)
TARGET_LINK_LIBRARIES(JDQ_Curve_Fusion zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin JDQ_Curve_Fusion)
