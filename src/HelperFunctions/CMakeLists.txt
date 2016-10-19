CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(HelperFunctions)


SET(HelperFunctions_INC
  inc/MylibValues.h
  inc/VectorFunctions.h
)

SET(HelperFunctions_SRC
  src/MylibValues.cpp
  src/VectorFunctions.cpp
)

SOURCE_GROUP("Header Files" FILES ${HelperFunctions_INC})
SOURCE_GROUP("Source Files" FILES ${HelperFunctions_SRC})

INCLUDE_DIRECTORIES(
  ${CMakeProject_HELPERFUNCTIONS_INC}
  ${CMakeProject_MYLIB_INC})

ADD_LIBRARY(HelperFunctions STATIC 
  ${HelperFunctions_INC}
  ${HelperFunctions_SRC}
)

TARGET_LINK_LIBRARIES(HelperFunctions
  ${CMakeProject_MYLIB}
)


