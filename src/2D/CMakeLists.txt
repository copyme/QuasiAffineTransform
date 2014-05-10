##########Check if libMagick++ exists
FIND_PACKAGE(Magick REQUIRED)

##tclap
INCLUDE_DIRECTORIES(${AQA_Transform_SOURCE_DIR}/src/tclap)

IF(MAGICK++_FOUND)
  INCLUDE_DIRECTORIES(${MAGICK++_INCLUDE_DIR})
  SET(optionalLibs ${optionalLibs} ${MAGICK++_LIBRARIES})
ENDIF(MAGICK++_FOUND)

SET(SRCs image matrix2x2 paving qat qat-2D vector2d ../common/mathematic)

ADD_EXECUTABLE(qat-2D ${SRCs})

TARGET_LINK_LIBRARIES(qat-2D ${optionalLibs})
