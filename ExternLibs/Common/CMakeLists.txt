#
# Source files
#
# Here is where you can add the name of your local Common classes.
#


PROJECT(ZXHCOMMON)

INCLUDE_DIRECTORIES(${ZXHCOMMON_SOURCE_DIR})

FILE(GLOB ZXHCommon_files "*.*")

SET (ZXHCommon_SRCS
${ZXHCommon_files}
)
 
# Create the zxhCommon C++ library
ADD_LIBRARY (zxhCommon ${ZXHCommon_SRCS})
 
TARGET_LINK_LIBRARIES(zxhCommon)
 
INSTALL_TARGETS(/lib/zxh zxhCommon)
INSTALL_FILES(/include/zxh .h ${ZXHCommon_SRCS})
