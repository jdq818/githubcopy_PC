#
# Source files
#
# Here is where you can add the name of your local ImageModelings classes.
#

PROJECT(zxhImageModelings)
INCLUDE_DIRECTORIES(${zxhImageModelings_SOURCE_DIR})

FILE(GLOB zxhImageModelings_files "*.h" "*.cpp")

SET (zxhImageModelings_SRCS
${zxhImageModelings_files}
)
 
# Create the zxhImageModelings C++ library
ADD_LIBRARY (zxhImageModelings ${zxhImageModelings_SRCS})
 

TARGET_LINK_LIBRARIES(zxhImageModelings zxhCommon zxhImages)

 
INSTALL_TARGETS(/lib  zxhImageModelings)
INSTALL_FILES(/include  .h ${zxhImageModelings_SRCS})
