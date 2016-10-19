CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(FlatFieldCorrection)

SET(FlatFieldCorrection_INC
)

SET(FlatFieldCorrection_SRC
  src/main.cpp
)

SOURCE_GROUP("Header Files" FILES ${FlatFieldCorrection_INC})
SOURCE_GROUP("Source Files" FILES ${FlatFieldCorrection_SRC})

INCLUDE_DIRECTORIES(
  ${CMakeProject_FLATFIELDCORRECTION_INC}
  ${CMakeProject_MYLIB_INC}
	${CMakeProject_HELPERFUNCTIONS_INC}
)

ADD_EXECUTABLE(FlatFieldCorrection 
	${FlatFieldCorrection_INC}
  ${FlatFieldCorrection_SRC})

TARGET_LINK_LIBRARIES(FlatFieldCorrection
  ${CMakeProject_MYLIB}
	${CMakeProject_HELPERFUNCTIONS}
)
