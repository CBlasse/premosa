cmake_minimum_required (VERSION 2.8)
project (mylib)
  
#
# Setup output directories for visual studio
#
if(MSVC OR XCODE)
  set(MYLIB_USE_GENERATOR_SPECIFIC_OUTPUT_DIR FALSE CACHE BOOL "Put built files into a directory named for the generator.  Helps support side-by-side 32-bit and 64-bit builds.")
  #message("!! ${CMAKE_CURRENT_SOURCE_DIR}")
  if(${MYLIB_USE_GENERATOR_SPECIFIC_OUTPUT_DIR})
    set(dir ${CMAKE_CURRENT_SOURCE_DIR}/${CMAKE_GENERATOR})
    set(EXECUTABLE_OUTPUT_PATH ${dir} CACHE PATH "Build directory" FORCE)
    set(LIBRARY_OUTPUT_PATH ${dir} CACHE PATH "Build directory" FORCE)
  endif()
endif()

#
# Useful tweaks
#
set(CMAKE_VERBOSE_MAKEFILE OFF)

if(NOT MSVC)
  add_definitions(-Wall -Wextra)
  # Define needed to get stat64 struct defined on redhat linux
  add_definitions(-D_LARGEFILE64_SOURCE)
  # On Linux AMD 64 systems, -fPIC is required to link static libraries to shared libraries
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC"  )
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC"  )
  set(CMAKE_LD_FLAGS "${CMAKE_LD_FLAGS} -fPIC"  )
endif()

#
# Configure
#

# check
#set(CMAKE_THREAD_PREFER_PTHREAD)
find_package(Threads)

#
# p-file transform
#
macro(show V)
message("${V} is ${${V}}")
endmacro()


add_executable(generator generator.c)
get_property(GENERATOR TARGET generator PROPERTY LOCATION )
macro(MYLIB_PFILE_TRANSFORM PFILES CFILES_OUT)
  file(GLOB ${PFILES} RELATIVE "${PROJECT_SOURCE_DIR}" *.p)
  foreach(f ${${PFILES}})
    string(REGEX REPLACE \\.p$ .c cfile ${f})
#message(":: ${f} --> ${cfile}")
    get_filename_component(src ${f} NAME)
    set(${CFILES_OUT} ${${CFILES_OUT}} ${PROJECT_BINARY_DIR}/${cfile})
    add_custom_command(
      OUTPUT ${PROJECT_BINARY_DIR}/${cfile}
      COMMAND ${GENERATOR} <${PROJECT_SOURCE_DIR}/${f} >${PROJECT_BINARY_DIR}/${cfile} 
      DEPENDS ${src} generator
      )
  endforeach(f)
  set_source_files_properties(${${CFILES_OUT}} PROPERTIES GENERATED 1)
endmacro()

#
# includes
#
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/MY_TIFF")
include_directories("${PROJECT_SOURCE_DIR}/MY_FFT")

#
# submodules
#
add_subdirectory (MY_TIFF)
add_subdirectory (MY_FFT)
set_source_files_properties(${MYTIFF_CFILES} PROPERTIES GENERATED 1)
set_source_files_properties(${MYFFT_GENERATED_CFILES}  PROPERTIES GENERATED 1)

#
# The Real Deal
#
file(GLOB MYLIB_P_SRCS *.p)
file(GLOB MYLIB_HDRS *.h)
MYLIB_PFILE_TRANSFORM(MYLIB_P_SRCS MYLIB_PtoC_SRCS)
file(GLOB MYLIB_C_SRCS *.c)
file(GLOB MYLIB_AWK_SRCS *.awk)

# Preserve default Makefile behavior of creating static library
set(DEFAULT_SHARED FALSE)
set(MYLIB_BUILD_SHARED_LIBRARY ${DEFAULT_SHARED} CACHE BOOL "Build shared rather than static library")
set(MYLIB_SRCS
  ${MYLIB_C_SRCS}
  ${MYLIB_P_SRCS}
  ${MYLIB_PtoC_SRCS}
  ${MYLIB_HDRS}
  ${MYLIB_AWK_SRCS}
)
if(MYLIB_BUILD_SHARED_LIBRARY)
  add_library(mylib MODULE ${MYLIB_SRCS})
else()
  add_library(mylib STATIC ${MYLIB_SRCS})
endif()

target_link_libraries(mylib mytiff)
target_link_libraries(mylib myfft)
target_link_libraries(mylib ${CMAKE_THREAD_LIBS_INIT})

#
# Source groups for IDE's that care (MSVC and XCode)
#

source_group("Header Files"          FILES ${MYLIB_HDRS})
source_group("Header Files\\MY_FFT"  FILES ${MYFFT_HDRS})
source_group("Header Files\\MY_TIFF" FILES ${MYTIFF_HFILES})
    
source_group("Source Files"          FILES ${MYLIB_C_SRCS} ${MYLIB_P_SRCS} ${MYLIB_AWK_SRCS})
source_group("Source Files\\MY_FFT"  FILES ${MYFFT_SRCS})
source_group("Source Files\\MY_TIFF" FILES ${MYTIFF_PFILES} ${MYTIFF_AWKFILES})
    
source_group("Source Files\\Generated"          FILES ${MYLIB_PtoC_SRCS})
source_group("Source Files\\Generated\\MY_FFT"  FILES ${MYFFT_GENERATED_CFILES})
source_group("Source Files\\Generated\\MY_TIFF" FILES ${MYTIFF_CFILES})
