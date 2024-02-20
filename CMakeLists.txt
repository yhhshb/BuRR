cmake_minimum_required(VERSION 3.11)
project(BuRR)

include(CMakePrintHelpers)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif ()

MESSAGE(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if (UNIX AND (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64"))
  MESSAGE(STATUS "Compiling with flags: -march=native -mavx -mbmi2 -msse4.2")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native") # If working on x86_64
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx") # for xxhash
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2 -msse4.2") # for hardware popcount and pdep
endif()

if (UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
  endif()
  if (BURR_USE_SANITIZERS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  endif()
endif()

include(FetchContent)
FetchContent_Declare(
  argparse
  GIT_REPOSITORY https://github.com/yhhshb/argparse.git
)
FetchContent_MakeAvailable(argparse)

include_directories(${PROJECT_SOURCE_DIR})

set(TLX_SOURCES
  bundled/tlx/src/cmdline_parser.cpp
  bundled/tlx/src/core.cpp
  bundled/tlx/src/parse_si_iec_units.cpp
)

set(BURR_SOURCES
  src/burr_config.cpp
  src/sorter.cpp
)

add_library(BuRR ${BURR_SOURCES})

if (MAKE_TESTS)
  add_executable(ribbon targets/main.cpp ${BURR_SOURCES} ${TLX_SOURCES})
  target_link_libraries(ribbon argparse rt)
endif()