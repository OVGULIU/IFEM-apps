project(LR-mesh)

cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_language(CXX)

find_package(LRSpline REQUIRED)

include_directories(${LRSpline_INCLUDE_DIRS})

set(APPS LRcrack LRcrack_lrline LRcrack_lrpoint
         LRcrack_shear LR_sneddon)

foreach(app ${APPS})
  add_executable(${app} ${app}.cpp)
  target_link_libraries(${app} ${LRSpline_LIBRARIES})
endforeach()
