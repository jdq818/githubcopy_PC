#
# Add the executable
#

ADD_EXECUTABLE(jdqvtk2txt jdqvtk2txt.cpp)
TARGET_LINK_LIBRARIES(jdqvtk2txt zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin jdqvtk2txt)
