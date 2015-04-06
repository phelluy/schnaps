set(CMAKE_LEGACY_CYGWIN_WIN32 0) # Remove when CMake >= 2.8.4 is required
cmake_minimum_required(VERSION 2.8)

project(schnaps)

#set(CMAKE_C_COMPILER "tau_cc.sh")
#set(CMAKE_C_FLAGS "-optRevert -optTauSelectFile=select.tau -std=c99 -fopenmp -Ofast")
#set(CMAKE_C_COMPILER "/opt/local/bin/gcc-mp-4.8")
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#set(CMAKE_C_FLAGS "-std=c99 -D_REENTRANT -O3 -g -D_POSIX_C_SOURCE=200112L")
#set(CMAKE_C_FLAGS "-std=gnu99 -D_REENTRANT -O3 -g")
set(CMAKE_C_FLAGS "-std=c99 -D_REENTRANT -O3 -g")

OPTION(USE_OPENCL
  "Compile with the OpenCL library"
  ON)

IF(USE_OPENCL)# set with: cmake . -DUSE_OPENCL:BOOL=ON
  MESSAGE("Compiling with OpenCL Library")
  find_package (OpenCL)
  if(OPENCL_FOUND)
    MESSAGE("OpenCL found")
  else()
    MESSAGE("OpenCL not found")
  endif(OPENCL_FOUND)
ENDIF(USE_OPENCL)

if(OPENCL_FOUND)
  INCLUDE_DIRECTORIES(${OPENCL_INCLUDE_DIR})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_OPENCL")
else()
  INCLUDE_DIRECTORIES(".")
endif(OPENCL_FOUND)

if(CMAKE_COMPILER_IS_GNUCC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror -Wno-deprecated-declarations")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
endif()

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(flann)

include_directories(src)

enable_testing()

list(APPEND
  source_files
  src/schnaps.h
  src/field.c
  src/field.h
  src/geometry.c
  src/geometry.h
  src/global.c
  src/global.h
  src/model.c
  src/model.h
  src/interpolation.c
  src/interpolation.h
  src/macromesh.c
  src/macromesh.h
  src/h20phi.h
  src/maxwell.h
  src/maxwell.c
  )

if (OPENCL_FOUND)
  list(APPEND
    source_files
    src/field_cl.c
    src/field_cl.h
    src/clutils.c
    src/clutils.h
    src/clinfo.c
    src/clinfo.h
    )
endif (OPENCL_FOUND)

#schnaps library
add_library (_schnaps ${source_files})
set(LIBS "")

if (OPENCL_FOUND)
list(APPEND LIBS _schnaps pthread ${OPENCL_LIBRARIES} m)
else()
list(APPEND LIBS _schnaps pthread m)
endif (OPENCL_FOUND)

# nearest neighbour fast search library
list(APPEND LIBS flann)


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

add_test(testmaxwell2d testmaxwell2d)
add_executable(
  testmaxwell2d
  test/test.h
  test/testmaxwell2d.c
  )
target_link_libraries(testmaxwell2d ${LIBS})
add_test(testmaxwell2d testmaxwell2d)

##########  OpenCL unit tests  ##########
if(OPENCL_FOUND)
  add_test(testclinfo testclinfo)
  add_executable(
    testclinfo
    test/test.h
    test/testclinfo.c
    )
  target_link_libraries(testclinfo ${LIBS})

  add_test(testmanyv testmanyv)
  add_executable(
    testmanyv
    test/test.h
    test/testmanyv.c
    )
  target_link_libraries(testmanyv ${LIBS})
  add_test(testmanyv testmanyv)


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

  add_test(testfieldrk4_cl testfieldrk4_cl)
  add_executable(
    testfieldrk4_cl
    test/test.h
    test/testfieldrk4_cl.c
    )
  target_link_libraries(testfieldrk4_cl ${LIBS})

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

endif(OPENCL_FOUND)
