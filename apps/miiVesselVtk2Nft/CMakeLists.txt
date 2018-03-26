#
# Add the executable
#

ADD_EXECUTABLE(miiVesselVtk2Nft miiVesselVtk2Nft.cpp)
TARGET_LINK_LIBRARIES(miiVesselVtk2Nft zxhImages zxhCommon znz zlib niftiio vtkRendering vtkIO)

INSTALL_TARGETS(/bin miiVesselVtk2Nft)
