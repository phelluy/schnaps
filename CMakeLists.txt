cmake_minimum_required(VERSION 2.8)

project(schnaps)

#set(CMAKE_C_COMPILER "tau_cc.sh")
#set(CMAKE_C_FLAGS "-optRevert -optTauSelectFile=select.tau -std=c99 -fopenmp -Ofast")
#set(CMAKE_C_COMPILER "/opt/local/bin/gcc-mp-4.8")
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
FIND_PACKAGE(OpenCL REQUIRED)
INCLUDE_DIRECTORIES(${OPENCL_INCLUDE_DIR})
#set(CMAKE_C_FLAGS "-std=c99 -D_REENTRANT -O3 -g -D_POSIX_C_SOURCE=200112L")
set(CMAKE_C_FLAGS "-std=gnu99 -D_REENTRANT -O3 -g")
#set(CMAKE_C_FLAGS "-std=c99 -D_REENTRANT -O3 -g")

if(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror -Wno-deprecated-declarations")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
endif()

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories(src)

enable_testing()

list(APPEND
source_files
src/clinfo.h
src/field.h
src/schnaps.h
src/geometry.h
src/global.c
src/global.h
src/h20phi.h
src/interpolation.h
src/macromesh.h
src/model.h
src/clinfo.c
src/field.c
src/geometry.c
src/interpolation.c
src/macromesh.c
src/model.c
src/clutils.c
src/clutils.h
)

#schnaps library
add_library (_schnaps ${source_files})
set(LIBS "")
list(APPEND LIBS _schnaps pthread ${OPENCL_LIBRARIES} rt m)

# main program schnaps
add_executable(
  schnaps
  src/schnaps.c
)
target_link_libraries(schnaps ${LIBS})

# list of tests
add_executable(
  testmacromesh
  test/test.h
  test/testmacromesh.c
)
target_link_libraries(testmacromesh ${LIBS})
add_test(testmacromesh testmacromesh)

add_executable(
  testgeometry
  test/test.h
  test/testgeometry.c
)
target_link_libraries(testgeometry ${LIBS})
add_test(testgeometry testgeometry)

add_executable(
  testinterpolation
  test/test.h
  test/testinterpolation.c
)
target_link_libraries(testinterpolation ${LIBS})
add_test(testinterpolation testinterpolation)

add_executable(
  testmodel
  test/test.h
  test/testmodel.c
)
target_link_libraries(testmodel ${LIBS})

add_test(testfield testfield)
add_executable(
  testfield
  test/test.h
  test/testfield.c
)
target_link_libraries(testfield ${LIBS})
add_test(testfield testfield)

add_test(testfielddg testfielddg)
add_executable(
  testfielddg
  test/test.h
  test/testfielddg.c
)
target_link_libraries(testfielddg ${LIBS})
add_test(testfielddg testfielddg)

add_test(testfieldrk2 testfieldrk2)
add_executable(
  testfieldrk2
  test/test.h
  test/testfieldrk2.c
)
target_link_libraries(testfieldrk2 ${LIBS})
add_test(testfieldrk2 testfieldrk2)

add_test(testfieldrk4 testfieldrk4)
add_executable(
  testfieldrk4
  test/test.h
  test/testfieldrk4.c
)
target_link_libraries(testfieldrk4 ${LIBS})
add_test(testfieldrk4 testfieldrk4)

add_test(test2dmeshdetection test2dmeshdetection)
add_executable(
  test2dmeshdetection
  test/test.h
  test/test2dmeshdetection.c
)
target_link_libraries(test2dmeshdetection ${LIBS})
add_test(test2dmeshdetection test2dmeshdetection)

add_test(testfieldrk2-2d testfieldrk2-2d)
add_executable(
  testfieldrk2-2d
  test/test.h
  test/testfieldrk2-2d.c
)
target_link_libraries(testfieldrk2-2d ${LIBS})
add_test(testfieldrk2-2d testfieldrk2-2d)

add_test(testsubcelldgvol testsubcelldgvol)
add_executable(
  testsubcelldgvol
  test/test.h
  test/testsubcelldgvol.c
)
target_link_libraries(testsubcelldgvol ${LIBS})
add_test(testsubcelldgvol testsubcelldgvol)

add_test(testfieldrk2-2d-subcell testfieldrk2-2d-subcell)

add_executable(
  testfieldrk2-2d-subcell
  test/test.h
  test/testfieldrk2-2d-subcell.c
)
target_link_libraries(testfieldrk2-2d-subcell ${LIBS})
add_test(testfieldrk2-2d-subcell testfieldrk2-2d-subcell)

add_test(testmeq2 testmeq2)
add_executable(
  testmeq2
  test/test.h
  test/testmeq2.c
)
target_link_libraries(testmeq2 ${LIBS})
add_test(testmeq2 testmeq2)

add_test(testmanyv testmanyv)
add_executable(
  testmanyv
  test/test.h
  test/testmanyv.c
)
target_link_libraries(testmanyv ${LIBS})
add_test(testmanyv testmanyv)

add_test(testclinfo testclinfo)
add_executable(
  testclinfo
  test/test.h
  test/testclinfo.c
)
target_link_libraries(testclinfo ${LIBS})

add_test(testkernel testkernel)
add_executable(
  testkernel
  test/test.h
  test/testkernel.c
)
target_link_libraries(testkernel ${LIBS})

add_test(testkernelvolume testkernelvolume)
add_executable(
  testkernelvolume
  test/test.h
  test/testkernelvolume.c
)
target_link_libraries(testkernelvolume ${LIBS})


add_test(testkernelinterface testkernelinterface)
add_executable(
  testkernelinterface
  test/test.h
  test/testkernelinterface.c
)
target_link_libraries(testkernelinterface ${LIBS})

add_test(testfieldrk2_cl testfieldrk2_cl)
add_executable(
  testfieldrk2_cl
  test/test.h
  test/testfieldrk2_cl.c
)
target_link_libraries(testfieldrk2_cl ${LIBS})

add_test(testdtfield_cl testdtfield_cl)
add_executable(
  testdtfield_cl
  test/test.h
  test/testdtfield_cl.c
)
target_link_libraries(testdtfield_cl ${LIBS})

add_test(testmacrocellinterface testmacrocellinterface)
add_executable(
  testmacrocellinterface
  test/test.h
  test/testmacrocellinterface.c
)
target_link_libraries(testmacrocellinterface ${LIBS})
