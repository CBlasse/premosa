CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(Projection)

SET(Projection_INC
)

SET(Projection_SRC
  src/main.cpp
)

SOURCE_GROUP("Header Files" FILES ${Projection_INC})
SOURCE_GROUP("Source Files" FILES ${Projection_SRC})

INCLUDE_DIRECTORIES(
  ${CMakeProject_PROJECTION_INC}
  ${CMakeProject_MYLIB_INC}
  ${CMakeProject_PROJECTIONLIBRARY_INC}
	${CMakeProject_HELPERFUNCTIONS_INC}
)

ADD_EXECUTABLE(Projection
  ${Projection_INC}
  ${Projection_SRC})

TARGET_LINK_LIBRARIES(Projection
  ${CMakeProject_MYLIB}
	${CMakeProject_HELPERFUNCTIONS}
	${CMakeProject_PROJECTIONLIBRARY}
)
