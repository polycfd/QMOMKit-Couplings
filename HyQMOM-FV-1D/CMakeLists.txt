cmake_minimum_required(VERSION 3.16)

project(hyqmomTransport1D LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# QMOMKIT_DIR is the QMOMKit install prefix, given on the command line
# (-DQMOMKIT_DIR=...) or in the environment.
find_package(QMOMKit 1.1 CONFIG REQUIRED
  HINTS
    ${QMOMKIT_DIR}
    $ENV{QMOMKIT_DIR}
  )

add_executable(hyqmomTransport1D
  hyqmomTransport1D.cc
  cases.cc
  exact_solution.cc
  fluxes.cc
  output.cc
  reference.cc
  settings.cc
  )
target_compile_features(hyqmomTransport1D PRIVATE cxx_std_17)
target_link_libraries(hyqmomTransport1D PRIVATE qmomkit::qmomkit)
