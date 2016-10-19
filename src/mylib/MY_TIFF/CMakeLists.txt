cmake_minimum_required(VERSION 2.8)
# OUTPUT
#
# MYTIFF_CFILES - a list of the generated files
# MYTIFF_HFILES - a list of the headers
#

file(GLOB AWK_SRCS    *.awk)
file(GLOB CORE_PFILES tiff.*.p)
file(GLOB CORE_HFILES tiff.*.h)
file(GLOB UTIL_PFILES utilities.p)
file(GLOB UTIL_HFILES utilities.h)

#
# Config
#

set(CMAKE_THREAD_PREFER_PTHREAD)
find_package(Threads)

#
# p-file transform
#
add_executable(mytiff_generator generator.c)
get_property(GENERATOR TARGET generator PROPERTY LOCATION )      
macro(MYTIFF_PFILE_TRANSFORM PFILES CFILES_OUT)
  foreach(f ${${PFILES}})
    get_filename_component(src ${f} NAME)
    string(REGEX REPLACE \\.p$ .c cfile ${src})
    #message(":: ${f} --> ${cfile}")
    #message("::      --> ${src}")
    set(${CFILES_OUT} ${${CFILES_OUT}} ${CMAKE_CURRENT_BINARY_DIR}/${cfile})
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${cfile}
      COMMAND ${GENERATOR} <${CMAKE_CURRENT_SOURCE_DIR}/${src} >${CMAKE_CURRENT_BINARY_DIR}/${cfile} 
      DEPENDS 
        ${src}
        mytiff_generator
      )
    #message("::      --> ${CMAKE_CURRENT_SOURCE_DIR}/${src}")
  endforeach(f)
  set_source_files_properties(${${CFILES_OUT}} PROPERTIES GENERATED 1)
endmacro()

##
# TARGETS
#
MYTIFF_PFILE_TRANSFORM(CORE_PFILES CORE_CFILES)
MYTIFF_PFILE_TRANSFORM(UTIL_PFILES UTIL_CFILES)

# foreach(f ${CORE_PFILES})
#   message("## ${f}")
# endforeach()
# message("###########")
# foreach(f ${CORE_CFILES})
#   message("## ${f}")
# endforeach()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

add_library(mytiff STATIC
    ${CORE_CFILES}
    ${CORE_PFILES}
    ${CORE_HFILES}
    ${AWK_SRCS}
    )

add_executable(tiffconvert 
    app.convert.c 
    ${UTIL_HFILES}
    ${UTIL_CFILES}
    ${UTIL_PFILES}
    )
target_link_libraries(tiffconvert mytiff)
target_link_libraries(tiffconvert ${CMAKE_THREAD_LIBS_INIT})

add_executable(tifftagger 
    app.tagger.c 
    ${UTIL_HFILES}
    ${UTIL_CFILES}
    ${UTIL_PFILES}
    )
target_link_libraries(tifftagger mytiff)
target_link_libraries(tifftagger ${CMAKE_THREAD_LIBS_INIT})

add_executable(tiffRGB 
    app.color.c 
    ${UTIL_HFILES}
    ${UTIL_CFILES}
    ${UTIL_PFILES}
    )
target_link_libraries(tiffRGB mytiff) 
target_link_libraries(tiffRGB ${CMAKE_THREAD_LIBS_INIT})

add_executable(tiffshow 
    app.show.c 
    ${UTIL_HFILES}
    ${UTIL_CFILES}
    ${UTIL_PFILES}
    )
target_link_libraries(tiffshow mytiff)
target_link_libraries(tiffshow ${CMAKE_THREAD_LIBS_INIT})

##
# Source groups (for IDEs that care)
#

source_group("Header Files" FILES
    ${CORE_HFILES}
    ${UTIL_HFILES}
    )
source_group("Source Files" FILES
    ${AWK_SRCS}
    app.show.c 
    app.color.c
    app.tagger.c
    app.convert.c
    ${UTIL_PFILES}
    ${CORE_PFILES}
    )    
source_group("Source Files\\Generated" FILES
    ${UTIL_CFILES}
    ${CORE_CFILES}
    )
##
# Output variables
#

set(MYTIFF_PFILES ${CORE_PFILES} PARENT_SCOPE)
set(MYTIFF_CFILES ${CORE_CFILES} PARENT_SCOPE)
set(MYTIFF_HFILES ${CORE_HFILES} PARENT_SCOPE)
set(MYTIFF_AWKFILES ${AWK_SRCS} PARENT_SCOPE)
