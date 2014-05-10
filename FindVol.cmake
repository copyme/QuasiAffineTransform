#-*-cmake-*-
#
# Test for LibVol libraries
#
# Once loaded this will define
#  VOL_FOUND        - system has libvol
#  VOL_INCLUDE_DIR  - include directory
#  VOL_LIBRARY_DIR  - library directory
#  VOL_LIBRARIES    - libraries you need to link to
#

SET(VOL_FOUND   "NO" )

FIND_PATH( VOL_INCLUDE_DIR vol.h
  "$ENV{VOL_LOCATION}"
  "$ENV{VOL_LOCATION}/include"
  "$ENV{VOL_HOME}/include"
  /usr/include/
  /usr/local/include/
  )


FIND_LIBRARY(Vol vol
  PATHS 
  "$ENV{VOL_LOCATION}/"
  "$ENV{VOL_LOCATION}/lib"
  "$ENV{VOL_HOME}/lib"
  DOC   "libvol library"
)



SET(VOL_LIBRARIES ${Vol} )


IF (VOL_INCLUDE_DIR)
  IF(VOL_LIBRARIES)
    SET(VOL_FOUND "YES")
    GET_FILENAME_COMPONENT(VOL_LIBRARY_DIR ${Vol}   PATH)
  ENDIF(VOL_LIBRARIES)
ENDIF(VOL_INCLUDE_DIR)



IF(NOT VOL_FOUND)
  # make FIND_PACKAGE friendly
  IF(NOT Vol_FIND_QUIETLY)
    IF(Vol_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR
              "libvol required, please specify it's location with VOL_HOME, VOL_LOCATION")
    ELSE(Vol_FIND_REQUIRED)
      MESSAGE(STATUS "libvol was not found.")
    ENDIF(Vol_FIND_REQUIRED)
  ENDIF(NOT Vol_FIND_QUIETLY)
ENDIF(NOT VOL_FOUND)


#####

