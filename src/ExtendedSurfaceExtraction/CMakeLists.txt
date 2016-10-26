CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(ExtendedSurfaceExtraction)

SET(ExtendedSurfaceExtraction_INC
  inc/Analysis.h
  inc/Develop.h
  inc/Filter.h
  inc/HelperFunctions.h
  inc/ImageClass.h
  inc/Interpolation.h
  inc/Refinement.h
  inc/ZPlaneSelection.h
  inc/ZProjection.h)

SET(ExtendedSurfaceExtraction_SRC
  src/main.cpp
  src/Analysis.cpp
  src/Filter.cpp
  src/HelperFunctions.cpp
  src/ImageClass.cpp
  src/Interpolation.cpp
  src/Refinement.cpp
  src/ZPlaneSelection.cpp
  src/ZProjection.cpp)

SOURCE_GROUP("Header Files" FILES ${ExtendedSurfaceExtraction_INC})
SOURCE_GROUP("Source Files" FILES ${ExtendedSurfaceExtraction_SRC})

INCLUDE_DIRECTORIES(
  ${CMakeProject_EXTENDEDSURFACEEXTRACTION_INC}
  ${CMakeProject_MYLIB_INC}
	${CMakeProject_HELPERFUNCTIONS_INC}
  ${CMakeProject_PROJECTIONLIBRARY_INC}
)

ADD_EXECUTABLE(ExtendedSurfaceExtraction
  ${ExtendedSurfaceExtraction_INC}
  ${ExtendedSurfaceExtraction_SRC})

TARGET_LINK_LIBRARIES(ExtendedSurfaceExtraction
  ${CMakeProject_MYLIB}
	${CMakeProject_HELPERFUNCTIONS}
  ${CMakeProject_PROJECTIONLIBRARY}
)
