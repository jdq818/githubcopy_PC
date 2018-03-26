
PROJECT(ZXHZLIB)
INCLUDE_REGULAR_EXPRESSION("^(deflate|inf|trees|zconf|zlib|zutil).*$")

INCLUDE_DIRECTORIES(${ZXHZLIB_SOURCE_DIR})
INCLUDE_DIRECTORIES(${ZXHZLIB_BINARY_DIR})

  
ADD_DEFINITIONS(-DHAVE_ZLIB)
INCLUDE_DIRECTORIES(zlib)
INCLUDE_DIRECTORIES(${ZXHZLIB_SOURCE_DIR})
INCLUDE_DIRECTORIES(${ZXHZLIB_BINARY_DIR})
SET(NIFTI_ZLIB_LIBRARIES zlib ${NIFTI_ZLIB_LIBRARIES})

# source files for zlib
FILE(GLOB ZXHZLIB_SRCS_Hfiles "*.h")
SET(ZXHZLIB_SRCS
adler32.c   gzio.c      inftrees.c  uncompr.c
compress.c  infblock.c  infutil.c   zutil.c
crc32.c     infcodes.c 
deflate.c   inffast.c  
inflate.c   trees.c 
${ZXHZLIB_SRCS_Hfiles}
)

# for windows add the .def and .rc files to the source list   // remove for building x64 in windows
#IF(WIN32)
#  IF(BUILD_SHARED_LIBS)
#    SET(ZLIB_DLL 1)
#  ENDIF(BUILD_SHARED_LIBS)
#  IF(NOT UNIX)
#    IF(NOT BORLAND)
#      IF(NOT MINGW)
#		SET(ZXHZLIB_SRCS ${ZXHZLIB_SRCS} zlib.def zlib.rc  )
#      ENDIF(NOT MINGW)
#    ENDIF(NOT BORLAND) 
#  ENDIF(NOT UNIX)   
#ENDIF(WIN32)

CONFIGURE_FILE(${ZXHZLIB_SOURCE_DIR}/.NoDartCoverage
			   ${ZXHZLIB_BINARY_DIR}/.NoDartCoverage)
CONFIGURE_FILE(${ZXHZLIB_SOURCE_DIR}/zlibDllConfig.h.in
			   ${ZXHZLIB_BINARY_DIR}/zlibDllConfig.h)


ADD_LIBRARY(zlib ${ZXHZLIB_SRCS})
target_link_libraries(zlib)

INSTALL_TARGETS(/lib/zxh zlib)
INSTALL_FILES(/include/zxh .h zlib zconf zlibDllConfig)

