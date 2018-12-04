#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.6)
#Name your project here
project(go)

#Sends the -std=c99 -lm flag to the gcc compiler
set(CMAKE_C_FLAGS "-std=gnu99 -g -O -Wall")

list(APPEND
  source_files
  dg_implicit.h
  lobatto.h
  skyline.h
  lobatto.c
  skyline.c
  )

#This tells CMake to fib.c and name it fibonacci
add_executable(go
  dg_implicit.c
  dg_implicit.h
  lobatto.h
  skyline.h
  lobatto.c
  skyline.c
  )

target_link_libraries (go m)

# cmake . -G "CodeLite - Unix Makefiles"
