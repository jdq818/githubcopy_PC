#
# Add the executable
#

ADD_EXECUTABLE(miiMinPathOrgApp miiMinPathOrgApp.cpp)
TARGET_LINK_LIBRARIES(miiMinPathOrgApp zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio )

INSTALL_TARGETS(/bin miiMinPathOrgApp)
