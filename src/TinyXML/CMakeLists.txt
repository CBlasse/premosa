CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(tinycml)


SET(TinyXML_INC
  inc/tinystr.h
  inc/tinyxml.h
)

SET(TinyXML_SRC
  src/tinystr.cpp
  src/tinyxml.cpp
  src/tinyxmlerror.cpp
  src/tinyxmlparser.cpp
)

SOURCE_GROUP("Header Files" FILES ${TinyXML_INC})
SOURCE_GROUP("Source Files" FILES ${TinyXML_SRC})

INCLUDE_DIRECTORIES(
  ${CMakeProject_TINYXML_INC}
  ${CMakeProject_MYLIB_INC})

ADD_LIBRARY(TinyXML STATIC 
  ${TinyXML_INC}
  ${TinyXML_SRC}
)

TARGET_LINK_LIBRARIES(TinyXML
  ${CMakeProject_MYLIB}
)


