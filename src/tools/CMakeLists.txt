##########Check if libMagick++ exists
FIND_PACKAGE(Magick REQUIRED)

##tclap
INCLUDE_DIRECTORIES(${AQA_Transform_SOURCE_DIR}/src/tclap)

IF(MAGICK++_FOUND)
  INCLUDE_DIRECTORIES(${MAGICK++_INCLUDE_DIR})
  SET(optionalLibs ${optionalLibs} ${MAGICK++_LIBRARIES})
ENDIF(MAGICK++_FOUND)

SET(SRCs QATPhantom2D ../2D/paving ../2D/image ../2D/vector2d ../common/mathematic ../common/color)

ADD_EXECUTABLE(QATPhantom2D ${SRCs})

TARGET_LINK_LIBRARIES(QATPhantom2D ${optionalLibs})
