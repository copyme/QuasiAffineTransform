# The name of our project is "AQA_Transform".
# David Coeurjolly <david.coeurjolly@liris.cnrs.fr>  08/2006

PROJECT(AQA_Transform)
SET(PROJECT_VERSION 0.2)

cmake_minimum_required(VERSION 2.4)

SET(CPACK_PACKAGE_VERSION 0.2)
SET(CPACK_PACKAGE_VERSION_MAJOR 	 0)
SET(CPACK_PACKAGE_VERSION_MINOR 	 2)
INCLUDE(CPack)

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/")

add_subdirectory(src)

SET(optionalLibs ${optionalLibs} m)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})


##Generation of the documentation using doxygen
#
#  new "make doc" target
#
INCLUDE("${CMAKE_MODULE_PATH}/TargetDoxygenDoc.cmake" OPTIONAL)
