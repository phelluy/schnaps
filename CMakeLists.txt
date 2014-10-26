cmake_minimum_required(VERSION 2.8)

project(schnaps)

set(CMAKE_C_FLAGS "-std=c99 -fopenmp -Ofast")
include_directories(src)

enable_testing()

list(APPEND
source_files
src/field.h
src/schnaps.h
src/geometry.h
src/global.h
src/h20phi.h
src/interpolation.h
src/macromesh.h
src/model.h
src/field.c
src/geometry.c
src/interpolation.c
src/macromesh.c
src/model.c
)

#schnaps library
add_library (_schnaps ${source_files})

# main program schnaps
add_executable(
  schnaps
  src/schnaps.c
)
target_link_libraries(schnaps _schnaps m)

# list of tests
add_executable(
  testmacromesh
  test/test.h
  test/testmacromesh.c
)
target_link_libraries(testmacromesh _schnaps m)
add_test(testmacromesh testmacromesh)

add_executable(
  testgeometry
  test/test.h
  test/testgeometry.c
)
target_link_libraries(testgeometry _schnaps m)
add_test(testgeometry testgeometry)

add_executable(
  testinterpolation
  test/test.h
  test/testinterpolation.c
)
target_link_libraries(testinterpolation _schnaps m)
add_test(testinterpolation testinterpolation)

add_executable(
  testmodel
  test/test.h
  test/testmodel.c
)
target_link_libraries(testmodel _schnaps m)

add_test(testfield testfield)
add_executable(
  testfield
  test/test.h
  test/testfield.c
)
target_link_libraries(testfield _schnaps m)
add_test(testfield testfield)

add_test(testfielddg testfielddg)
add_executable(
  testfielddg
  test/test.h
  test/testfielddg.c
)
target_link_libraries(testfielddg _schnaps m)
add_test(testfielddg testfielddg)

add_test(testfieldrk2 testfieldrk2)
add_executable(
  testfieldrk2
  test/test.h
  test/testfieldrk2.c
)
target_link_libraries(testfieldrk2 _schnaps m)
add_test(testfieldrk2 testfieldrk2)

add_test(test2dmeshdetection test2dmeshdetection)
add_executable(
  test2dmeshdetection
  test/test.h
  test/test2dmeshdetection.c
)
target_link_libraries(test2dmeshdetection _schnaps m)
add_test(test2dmeshdetection test2dmeshdetection)

add_test(testfieldrk2-2d testfieldrk2-2d)
add_executable(
  testfieldrk2-2d
  test/test.h
  test/testfieldrk2-2d.c
)
target_link_libraries(testfieldrk2-2d _schnaps m)
add_test(testfieldrk2-2d testfieldrk2-2d)


add_test(testsubcelldgvol testsubcelldgvol)
add_executable(
  testsubcelldgvol
  test/test.h
  test/testsubcelldgvol.c
)
target_link_libraries(testsubcelldgvol _schnaps m)
add_test(testsubcelldgvol testsubcelldgvol)

add_test(testfieldrk2-2d-subcell testfieldrk2-2d-subcell)
add_executable(
  testfieldrk2-2d-subcell
  test/test.h
  test/testfieldrk2-2d-subcell.c
)
target_link_libraries(testfieldrk2-2d-subcell _schnaps m)
add_test(testfieldrk2-2d-subcell testfieldrk2-2d-subcell)
