#-*-cmake-*-
#
# Test for LibVol libraries
#
# Once loaded this will define
#  LONGVOL_FOUND        - system has libvol
#  LONGVOL_INCLUDE_DIR  - include directory
#  LONGVOL_LIBRARY_DIR  - library directory
#  LONGVOL_LIBRARIES    - libraries you need to link to
#

SET(LONGVOL_FOUND   "NO" )

FIND_PATH( LONGVOL_INCLUDE_DIR longvol.h
  "$ENV{LONGVOL_LOCATION}"
  "$ENV{LONGVOL_LOCATION}/include"
  "$ENV{LONGVOL_HOME}/include"
  /usr/include/
  /usr/local/include/
  )


FIND_LIBRARY(Longvol longvol
  PATHS 
  "$ENV{LONGVOL_LOCATION}/"
  "$ENV{LONGVOL_LOCATION}/lib"
  "$ENV{LONGVOL_HOME}/lib"
  DOC   "libvol library"
)



SET(LONGVOL_LIBRARIES ${Longvol} )


IF (LONGVOL_INCLUDE_DIR)
  IF(LONGVOL_LIBRARIES)
    SET(LONGVOL_FOUND "YES")
    GET_FILENAME_COMPONENT(LONGVOL_LIBRARY_DIR ${Longvol}   PATH)
  ENDIF(LONGVOL_LIBRARIES)
ENDIF(LONGVOL_INCLUDE_DIR)



IF(NOT LONGVOL_FOUND)
  # make FIND_PACKAGE friendly
  IF(NOT Vol_FIND_QUIETLY)
    IF(Vol_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR
              "libvol required, please specify it's location with LONGVOL_HOME, LONGVOL_LOCATION")
    ELSE(Vol_FIND_REQUIRED)
      MESSAGE(STATUS "liblongvol was not found.")
    ENDIF(Vol_FIND_REQUIRED)
  ENDIF(NOT Vol_FIND_QUIETLY)
ENDIF(NOT LONGVOL_FOUND)


#####

