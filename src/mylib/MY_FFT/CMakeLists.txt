cmake_minimum_required(VERSION 2.8)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..)

file(GLOB FFT_HDRS *.h)
file(GLOB FFTF fft.F.c) 
file(GLOB FFTD fft.D.c) 
file(GLOB FFTB myfft.c) 

source_group("Header Files" FILES ${FFT_HDRS})
source_group("Source Files" FILES ${FFTF} ${FFTD} ${FFTB})

set(MYFFT_SRCS ${FFTF} ${FFTD} ${FFTB} PARENT_SCOPE)
set(MYFFT_HDRS ${FFT_HDRS} PARENT_SCOPE)

add_library(myfft STATIC
    ${FFT_HDRS}
    ${FFTF}
    ${FFTD}
    ${FFTDB}
    )
