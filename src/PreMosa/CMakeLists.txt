CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(PreMosa)

SET(PreMosa_INC
	inc/FileHandling.h
	inc/ParseXML.h
	inc/Pipeline.h
  )

SET(PreMosa_SRC
  src/main.cpp
	src/FileHandling.cpp
	src/ParseXML.cpp
	src/Pipeline.cpp
  )

SOURCE_GROUP("Header Files" FILES ${PreMosa_INC})
SOURCE_GROUP("Source Files" FILES ${PreMosa_SRC})

INCLUDE_DIRECTORIES(
  ${CMakeProject_PREMOSA_INC}
  ${CMakeProject_MYLIB_INC}
	${CMakeProject_HELPERFUNCTIONS_INC}
  ${CMakeProject_PROJECTIONLIBRARY_INC}
)

ADD_EXECUTABLE(PreMosa
  ${PreMosa_INC}
  ${PreMosa_SRC})

TARGET_LINK_LIBRARIES(PreMosa
  ${CMakeProject_MYLIB}
	${CMakeProject_HELPERFUNCTIONS}
	${CMakeProject_PROJECTIONLIBRARY}
	Qt5::Core
	Qt5::Xml
)
