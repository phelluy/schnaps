cmake_minimum_required(VERSION 2.8)

project(schnaps)

set(CMAKE_C_FLAGS "-std=c99 -O -lm")

enable_testing()

list(APPEND
source_files
src/field.h
src/geometry.h
src/global.h
src/h20phi.h
src/interpolation.h
src/macromesh.h
src/model.h
src/test.h
src/field.c
src/geometry.c
src/interpolation.c
src/macromesh.c
src/model.c
)

# main program schnaps
add_executable(
  schnaps
  ${source_files}
  src/schnaps.c
)
target_link_libraries(schnaps m)

# list of tests
add_executable(
  testmacromesh
  ${source_files}
  src/testmacromesh.c
)
target_link_libraries(testmacromesh m)
add_test(testmacromesh testmacromesh)

add_executable(
  testgeometry
  ${source_files}
  src/testgeometry.c
)
target_link_libraries(testgeometry m)
add_test(testgeometry testgeometry)

add_executable(
  testinterpolation
  ${source_files}
  src/testinterpolation.c
)
target_link_libraries(testinterpolation m)
add_test(testinterpolation testinterpolation)

add_executable(
  testmodel
  ${source_files}
  src/testmodel.c
)
target_link_libraries(testmodel m)

add_test(testfield testfield)
add_executable(
  testfield
  ${source_files}
  src/testfield.c
)
target_link_libraries(testfield m)
add_test(testfield testfield)

add_test(testfielddg testfielddg)
add_executable(
  testfielddg
  ${source_files}
  src/testfielddg.c
)
target_link_libraries(testfielddg m)
add_test(testfielddg testfielddg)

add_test(testfieldrk2 testfieldrk2)
add_executable(
  testfieldrk2
  ${source_files}
  src/testfieldrk2.c
)
target_link_libraries(testfieldrk2 m)
add_test(testfieldrk2 testfieldrk2)
