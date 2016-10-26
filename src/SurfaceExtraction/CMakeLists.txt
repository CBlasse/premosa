CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(SurfaceExtraction)

SET(SurfaceExtraction_INC
)

SET(SurfaceExtraction_SRC
  src/main.cpp
)

SOURCE_GROUP("Header Files" FILES ${SurfaceExtraction_INC})
SOURCE_GROUP("Source Files" FILES ${SurfaceExtraction_SRC})

INCLUDE_DIRECTORIES(
  ${CMakeProject_SURFACEEXTRACTION_INC}
  ${CMakeProject_MYLIB_INC}
  ${CMakeProject_PROJECTIONLIBRARY_INC}
	${CMakeProject_HELPERFUNCTIONS_INC}
)

ADD_EXECUTABLE(SurfaceExtraction
  ${SurfaceExtraction_INC}
  ${SurfaceExtraction_SRC})

TARGET_LINK_LIBRARIES(SurfaceExtraction
  ${CMakeProject_MYLIB}
	${CMakeProject_HELPERFUNCTIONS}
	${CMakeProject_PROJECTIONLIBRARY}
)
