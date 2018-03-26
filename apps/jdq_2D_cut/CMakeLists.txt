#
# Add the executable
#

ADD_EXECUTABLE(jdq_2D_cut jdq_2D_cut.cpp)
TARGET_LINK_LIBRARIES(jdq_2D_cut zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO zxhImageModelings)

INSTALL_TARGETS(/bin jdq_2D_cut)
