cmake_minimum_required(VERSION 2.8)

project(schnaps)

set(CMAKE_C_FLAGS "-std=c99 -O")

file(
  GLOB_RECURSE
  source_files
  src/*.h
  src/*.c
)

add_executable(

  schnaps
  ${source_files}
)

target_link_libraries(schnaps m)
