#
# Add the executable
#

ADD_EXECUTABLE(jdqimageop_recipr jdqimageop_recipr.cpp)
TARGET_LINK_LIBRARIES(jdqimageop_recipr zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin jdqimageop_recipr)
