CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(ProjectionLibrary)

SET(ProjectionLibrary_INC
  inc/Develop.h
  inc/ExtendedProjectionAlgorithm.h
  inc/Filter.h
  inc/HeightMapAnalysis.h
  inc/Interpolation.h
	inc/Parameter.h
  inc/ProjectionAlgorithm.h
	inc/RasterizedImage.h
  inc/Refinement.h
  inc/ShowArray.h
  inc/ZPlaneSelection.h
  inc/ZProjection.h)

SET(ProjectionLibrary_SRC
  src/ExtendedProjectionAlgorithm.cpp
  src/Filter.cpp
  src/HeightMapAnalysis.cpp
  src/Interpolation.cpp
	src/Parameter.cpp
  src/ProjectionAlgorithm.cpp
	src/RasterizedImage.cpp
  src/Refinement.cpp
  src/ShowArray.cpp
  src/ZPlaneSelection.cpp
  src/ZProjection.cpp)

SOURCE_GROUP("Header Files" FILES ${ProjectionLibrary_INC})
SOURCE_GROUP("Source Files" FILES ${ProjectionLibrary_SRC})

INCLUDE_DIRECTORIES(
  ${CMakeProject_PROJECTIONLIBRARY_INC}
  ${CMakeProject_MYLIB_INC}
	${CMakeProject_HELPERFUNCTIONS_INC}
)

ADD_LIBRARY(ProjectionLibrary STATIC 
  ${ProjectionLibrary_INC}
  ${ProjectionLibrary_SRC}
)

TARGET_LINK_LIBRARIES(ProjectionLibrary
  ${CMakeProject_MYLIB}
	${CMakeProject_HELPERFUNCTIONS}
)
