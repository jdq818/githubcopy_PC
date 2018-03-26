#
# Add the executable
#

ADD_EXECUTABLE(miiFindMaxPath miiFindMaxPath.cpp)
TARGET_LINK_LIBRARIES(miiFindMaxPath zxhImages zxhCommon znz zlib niftiio vtkRendering vtkIO)

INSTALL_TARGETS(/bin miiFindMaxPath)
