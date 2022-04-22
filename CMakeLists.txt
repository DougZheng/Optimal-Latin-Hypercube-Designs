cmake_minimum_required(VERSION 3.5)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message("Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

project("LHD")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(LHD_algorithm
  "src/comm/design.cpp"
  "src/comm/design.h"
  "src/comm/utils.h"
  "src/algorithm/algorithm_manager.h"
  "src/algorithm/construction_algorithm.h"
  "src/algorithm/ESE.cpp"
  "src/algorithm/ESE.h"
  "src/algorithm/GA.cpp"
  "src/algorithm/GA.h"
  "src/algorithm/LaPSO.cpp"
  "src/algorithm/LaPSO.h"
  "src/algorithm/LSGA.cpp"
  "src/algorithm/LSGA.h"
  "src/algorithm/SA.cpp"
  "src/algorithm/SA.h"
  "src/algorithm/search_algorithm.h"
  "src/algorithm/Wang2018.cpp"
  "src/algorithm/Wang2018.h"
  "src/criterion/avgabscor_criterion.cpp"
  "src/criterion/avgabscor_criterion.h"
  "src/criterion/criterion.h"
  "src/criterion/maxabscor_criterion.cpp"
  "src/criterion/maxabscor_criterion.h"
  "src/criterion/multicriteria.cpp"
  "src/criterion/multicriteria.h"
  "src/criterion/phipl1_criterion.cpp"
  "src/criterion/phipl1_criterion.h"
  "src/criterion/phipl2_criterion.cpp"
  "src/criterion/phipl2_criterion.h"
)

include_directories(${PROJECT_SOURCE_DIR}/src)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../bin)

add_executable(LHD
  "src/LHD.cpp"
)

target_link_libraries(LHD LHD_algorithm)

add_executable(metric
  "src/metric.cpp"
)

add_executable(example
  "src/example.cpp"
)

target_link_libraries(example LHD_algorithm)