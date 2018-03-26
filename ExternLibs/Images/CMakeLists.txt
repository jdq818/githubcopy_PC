#
# Source files
#
# Here is where you can add the name of your local Transform classes.
#


PROJECT(zxhImageProcessings)
INCLUDE_DIRECTORIES(${zxhImageProcessings_SOURCE_DIR})


FILE(GLOB zxhImages_files "*.h" "*.cpp")

SET (zxhImages_SRCS
${zxhImages_files}
)
 
# Create the zxhImages C++ library
ADD_LIBRARY (zxhImages ${zxhImages_SRCS})


TARGET_LINK_LIBRARIES(zxhImages zxhCommon niftiio znz ${NIFTI_NIFTILIB_NAME} ${VTK_LIBS})

INSTALL_TARGETS(/lib  zxhImages)

INSTALL_FILES(/include  .h ${zxhImages_SRCS})
