##########Check if libVol exists
FIND_PACKAGE(Vol REQUIRED)

IF(VOL_FOUND)
  INCLUDE_DIRECTORIES(${VOL_INCLUDE_DIR})
  SET(optionalLibs ${optionalLibs} ${VOL_LIBRARIES})
ENDIF(VOL_FOUND)

##########Check if libLongVol exists
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/")
FIND_PACKAGE(Longvol REQUIRED)


##tclap
INCLUDE_DIRECTORIES(${AQA_Transform_SOURCE_DIR}/src/tclap)

IF(LONGVOL_FOUND)
  INCLUDE_DIRECTORIES(${LONGVOL_INCLUDE_DIR})
  SET(optionalLibs ${optionalLibs} ${LONGVOL_LIBRARIES})
ENDIF(LONGVOL_FOUND)

SET(SRCs qat-3D qat3d image3d paving3d vector3d ../common/mathematic)

ADD_EXECUTABLE(qat-3D ${SRCs})

TARGET_LINK_LIBRARIES(qat-3D ${optionalLibs})
