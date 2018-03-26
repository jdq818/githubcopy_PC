#
# Source files
#
# Here is where you can add the name of your local Transform classes.
#

FILE(GLOB ALL_FILES "*.*")

SET (ALL_SRCS
${ALL_FILES}
)
 
# Create the minimal path C++ library
ADD_LIBRARY (zxhcaeMinPathLib ${ALL_SRCS})


TARGET_LINK_LIBRARIES(zxhcaeMinPathLib zxhCommon niftiio znz ${NIFTI_NIFTILIB_NAME}  vtkRendering vtkIO)

INSTALL_TARGETS(/lib/zxh zxhcaeMinPathLib)

INSTALL_FILES(/include/zxh .h ${ALL_SRCS})
