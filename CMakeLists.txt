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

#schnaps library
add_library (libschnaps ${source_files})

# main program schnaps
add_executable(
  schnaps
  src/schnaps.c
)
target_link_libraries(schnaps libschnaps m)

# list of tests
add_executable(
  testmacromesh
  src/testmacromesh.c
)
target_link_libraries(testmacromesh libschnaps m)
add_test(testmacromesh testmacromesh)

add_executable(
  testgeometry
  src/testgeometry.c
)
target_link_libraries(testgeometry libschnaps m)
add_test(testgeometry testgeometry)

add_executable(
  testinterpolation
  src/testinterpolation.c
)
target_link_libraries(testinterpolation libschnaps m)
add_test(testinterpolation testinterpolation)

add_executable(
  testmodel
  src/testmodel.c
)
target_link_libraries(testmodel libschnaps m)

add_test(testfield testfield)
add_executable(
  testfield
  src/testfield.c
)
target_link_libraries(testfield libschnaps m)
add_test(testfield testfield)

add_test(testfielddg testfielddg)
add_executable(
  testfielddg
  src/testfielddg.c
)
target_link_libraries(testfielddg libschnaps m)
add_test(testfielddg testfielddg)

add_test(testfieldrk2 testfieldrk2)
add_executable(
  testfieldrk2
  src/testfieldrk2.c
)
target_link_libraries(testfieldrk2 libschnaps m)
add_test(testfieldrk2 testfieldrk2)
