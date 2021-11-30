cmake_minimum_required(VERSION 3.9)
set(PRGNAME MUSIC)
project(MUSIC)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -Wall -pedantic -DCMAKE_BUILD")
find_package(PkgConfig REQUIRED)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}")

option(MUSIC_ENABLE_SINGLE_PRECISION "Enable Single Precision Mode" OFF)

########################################################################################################################
# OpenMP
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

########################################################################################################################
# Pthreads
find_package(Threads REQUIRED)

########################################################################################################################
# FFTW
if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()
find_package(FFTW3 COMPONENTS SINGLE DOUBLE THREADS)

########################################################################################################################
# TIRPC, needed only for Tipsy format
find_package(TIRPC)

########################################################################################################################
# Add a custom command that produces version.cc, plus
# a dummy output that's not actually produced, in order
# to force version.cmake to always be re-run before the build
ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/version.cc
           ${CMAKE_CURRENT_BINARY_DIR}/_version.cc
    COMMAND ${CMAKE_COMMAND} -P
            ${CMAKE_CURRENT_SOURCE_DIR}/version.cmake)



########################################################################################################################
# GSL
find_package(GSL REQUIRED)

########################################################################################################################
# HDF5
find_package(HDF5)

########################################################################################################################
# INCLUDES
include_directories(${PROJECT_SOURCE_DIR}/src)

# SOURCES
# get all the *.cc files in the subfolders
file( GLOB SOURCES
  ${PROJECT_SOURCE_DIR}/src/*.cc
)
# add the auto generated version file
list (APPEND SOURCES "${CMAKE_CURRENT_BINARY_DIR}/version.cc")

# PLUGINS
# get all the *.cc files in the plugin subfolder
file( GLOB PLUGINS
  ${PROJECT_SOURCE_DIR}/src/plugins/*.cc
)

add_executable(${PRGNAME} ${SOURCES} ${PLUGINS})

set_target_properties(${PRGNAME} PROPERTIES CXX_STANDARD 11)

if(FFTW3_FOUND)
  target_compile_options(${PRGNAME} PRIVATE "-DFFTW3")

  if( MUSIC_ENABLE_SINGLE_PRECISION )
    target_compile_options(${PRGNAME} PRIVATE "-DSINGLE_PRECISION")
    if (FFTW3_SINGLE_THREADS_FOUND)
      target_link_libraries(${PRGNAME} ${FFTW3_SINGLE_THREADS_LIBRARY})
      target_compile_options(${PRGNAME} PRIVATE "-DUSE_FFTW_THREADS")
    elseif(FFTW3_SINGLE_SERIAL_FOUND)
      target_link_libraries(${PRGNAME} ${FFTW3_SINGLE_SERIAL_LIBRARY})
      message( WARNING "using serial version of FFTW3 -- this will most likely cause a very slow version of MUSIC. Rec: install FFTW3 with thread support")
    else()  
      message( FATAL "chose compilation in single precision, but FFTW3 not found for single precision")
    endif()
  else(MUSIC_ENABLE_SINGLE_PRECISION)
    if (FFTW3_DOUBLE_THREADS_FOUND)
      target_link_libraries(${PRGNAME} ${FFTW3_DOUBLE_THREADS_LIBRARY})
      target_compile_options(${PRGNAME} PRIVATE "-DUSE_FFTW_THREADS")
    elseif(FFTW3_DOUBLE_SERIAL_FOUND)
      target_link_libraries(${PRGNAME} ${FFTW3_DOUBLE_SERIAL_LIBRARY})
      message( WARNING "using serial version of FFTW3 -- this will most likely cause a very slow version of MUSIC. Rec: install FFTW3 with thread support")
    else()  
      message( FATAL "chose compilation in double precision, but FFTW3 not found for double precision")
    endif()
  endif(MUSIC_ENABLE_SINGLE_PRECISION)
endif(FFTW3_FOUND)

if(HDF5_FOUND)
  # target_link_libraries(${PRGNAME} ${HDF5_C_LIBRARY_DIRS})
  target_link_libraries(${PRGNAME} ${HDF5_LIBRARIES})
  target_include_directories(${PRGNAME} PRIVATE ${HDF5_INCLUDE_DIRS})
  target_compile_options(${PRGNAME} PRIVATE "-DHAVE_HDF5")
  target_compile_options(${PRGNAME} PRIVATE "-DH5_USE_16_API")
endif(HDF5_FOUND)

if(TIRPC_FOUND)
  target_link_libraries(${PRGNAME} ${TIRPC_LIBRARIES})
  target_include_directories(${PRGNAME} PRIVATE ${TIRPC_INCLUDE_DIRS})
  target_compile_options(${PRGNAME} PRIVATE "-DHAVE_TIRPC")
endif(TIRPC_FOUND)

target_link_libraries(${PRGNAME} ${FFTW3_LIBRARIES})
target_include_directories(${PRGNAME} PRIVATE ${FFTW3_INCLUDE_DIRS})

target_link_libraries(${PRGNAME} ${GSL_LIBRARIES})
target_include_directories(${PRGNAME} PRIVATE ${GSL_INCLUDE_DIR})

target_link_libraries(${PRGNAME} ${HDF5_LIBRARIES})
target_include_directories(${PRGNAME} PRIVATE ${HDF5_INCLUDE_DIR})
