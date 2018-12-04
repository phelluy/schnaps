cmake_minimum_required(VERSION 2.8.4)


# SET(CMAKE_VERBOSE_MAKEFILE ON)

# .................CMAKE BUILD TYPE.......................
IF(CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE}
    CACHE STRING "Choose the type of build, options are: Debug Release.")
ELSE()
  SET(CMAKE_BUILD_TYPE "Release"
    CACHE STRING "Choose the type of build, options are: Debug Release." FORCE)
ENDIF()
# .......................................................

# .......................PROJECT.........................
PROJECT(schnaps)
# .......................................................

# .......................................................
# Prevents developers to compile from the source directory
IF (${CMAKE_SOURCE_DIR} MATCHES ${CMAKE_BINARY_DIR})
  MESSAGE(SEND_ERROR)
  MESSAGE("-- ERROR: in-tree-build not allowed.")
  MESSAGE("-- TRY: rm -f CMakeCache.txt; mkdir build; cd build; cmake ..; make")
  RETURN()
ENDIF()
# .......................................................

#.....................CMAKE C FLAGS......................
SET(CMAKE_C_FLAGS " -std=gnu99 -D_SCHNAPS_DIR=${PROJECT_SOURCE_DIR} -O3 ")
 # SET(CMAKE_C_FLAGS " -g ")
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_DOUBLE_PRECISION -g")
# Adding warning flags:
IF(CMAKE_COMPILER_IS_GNUCC)
  #SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror -Wno-deprecated-declarations")
  #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
ENDIF(CMAKE_COMPILER_IS_GNUCC)

MESSAGE(STATUS "# Compiling options: ${CMAKE_C_FLAGS}")
# .......................................................

#..................CMAKE MODULE PATH.....................
# Add new directory to find the FIND<PACKAGE>.cmake files. CMAKE_SOURCE_DIR is
# the directory where the project sources are located. CMAKE_MODULE_PATH is 
# the list of directories in which to search for CMake modules.
SET (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH} )
# .......................................................

# ......................USE PYTHON3.......................
OPTION(USE_PYTHON "Activate Python Libs library" ON)
IF(USE_PYTHON) # set with: cmake . -DUSE_OPENCL:BOOL=ON
  MESSAGE("# Compiling with Python3 embedding Library...")
  FIND_PACKAGE (PythonLibs 3 QUIET)
  IF(PYTHONLIBS_FOUND)
    MESSAGE(STATUS "Python Libs found ----------")
    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIR})
    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
    GET_FILENAME_COMPONENT(PYTHON_LIBRARY_DIR ${PYTHON_LIBRARY} DIRECTORY)
    LINK_DIRECTORIES(${PYTHON3_LIBRARY_DIR})
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_PYTHONLIBS")
  ELSE()
    MESSAGE(STATUS "PYTHON NOT found ----------")
  ENDIF(PYTHONLIBS_FOUND)
ELSE()
  MESSAGE(STATUS "PYTHONLIBS_FOUND is not activated ----------")
ENDIF(USE_PYTHON)



# ......................USE OPENCL.......................
OPTION(USE_OPENCL "Activate OpenCL library" ON)

IF(USE_OPENCL) # set with: cmake . -DUSE_OPENCL:BOOL=ON
  MESSAGE("# Compiling with OpenCL Library...")
  FIND_PACKAGE (OpenCL)
  IF(OPENCL_FOUND)
    MESSAGE(STATUS "OpenCL found ----------")
    INCLUDE_DIRECTORIES(${OPENCL_INCLUDE_DIR})
    INCLUDE_DIRECTORIES(${OPENCL_INCLUDE_DIRS})
    GET_FILENAME_COMPONENT(OPENCL_LIBRARY_DIR ${OPENCL_LIBRARIES} DIRECTORY)
    LINK_DIRECTORIES(${OPENCL_LIBRARY_DIR})
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_OPENCL")
  ELSE()
    MESSAGE(STATUS "OpenCL NOT found ----------")
  ENDIF(OPENCL_FOUND)
ELSE()
  MESSAGE(STATUS "OpenCL is not activated ----------")
ENDIF(USE_OPENCL)
# .......................................................

# ......................USE STARPU.......................
# TODO: this should probably not be an option any more...
OPTION(USE_STARPU "Compile with the StarPU library" ON)

# NB: FindSTARPU.cmake does not properly deal with non-standard
# locations, so it may be necessary to just put the include directory
# in by hand.  For example,
# cmake ../ -DSTARPU_INCDIR=/home/<username>/starpu/include/starpu/1.2/
IF(USE_STARPU)
  MESSAGE("# Compiling with StarPU Library...")
  IF(STARPU_INCDIR)
    INCLUDE_DIRECTORIES(${STARPU_INCDIR})
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_STARPU")
  ELSE()
    FIND_PACKAGE (STARPU REQUIRED COMPONENTS HWLOC MPI)
    IF(STARPU_FOUND)
      MESSAGE(STATUS "StarPU found ----------")
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_STARPU")
      INCLUDE_DIRECTORIES(${STARPU_INCLUDE_DIRS})
      LINK_DIRECTORIES(${STARPU_LIBRARY_DIRS})
    ELSE()
      MESSAGE(STATUS "StarPU not found ----------")
    ENDIF(STARPU_FOUND)
  ENDIF(STARPU_INCDIR)
ENDIF(USE_STARPU)
# .......................................................

# ......................USE IGRAPH.......................
# TODO: this should probably not be an option any more...
SET(IGRAPH_ENABLED   ON CACHE BOOL   "Enable igraph use")

IF(IGRAPH_ENABLED)
  MESSAGE("# Compiling with igraph Library...")
  IF(IGRAPH_INCDIR)
    INCLUDE_DIRECTORIES(${IGRAPH_INCDIR})
  ELSE()
    FIND_PACKAGE (igraph)
    IF(HAVE_IGRAPH)
      MESSAGE(STATUS "igrah found ----------")
      INCLUDE_DIRECTORIES(${IGRAPH_INCLUDES})
    ELSE()
      MESSAGE(STATUS "igrah not found ----------")
    ENDIF(HAVE_IGRAPH)
  ENDIF(IGRAPH_INCDIR)
ENDIF(IGRAPH_ENABLED)
# .......................................................

# ................USE Garbage Collector..................
# TODO: Put default value to ON
SET(GC_ENABLED   OFF CACHE BOOL   "Enable Garbage Collector")

IF(GC_ENABLED)# set with: cmake . -DGC_ENABLED:BOOL=ON
  MESSAGE("# Compiling with Boehm GC...")
  FIND_PACKAGE (GC)
  IF(HAVE_GC)
    MESSAGE("GC found ----------")
    INCLUDE_DIRECTORIES(${BOEHM_GC_INCLUDE_DIR})
  ELSE()
    MESSAGE("GC not found ----------")
  ENDIF(HAVE_GC)
ENDIF(GC_ENABLED)
# .......................................................

# ......................SUITESPARSE......................
SET(SuiteSparse_VERBOSE OFF)
SET(SuiteSparse_USE_LAPACK_BLAS ON)
FIND_PACKAGE(suitesparse)
IF(SuiteSparse_FOUND)
  MESSAGE("# Compiling with SuiteSparse...")
  MESSAGE(STATUS "# umfpack libs ok:" ${SuiteSparse_INCLUDE_DIRS})
  LIST(APPEND INCLUDE_DIRECTORIES ${SuiteSparse_INCLUDE_DIRS})
  INCLUDE_DIRECTORIES(${SuiteSparse_INCLUDE_DIRS})
  LIST(APPEND LIBS ${SuiteSparse_LAPACK_BLAS_LIBRARIES})
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_UMFPACK")
ELSE()
  MESSAGE(STATUS "umfpack not found: deactivated ----------")
ENDIF()
# .......................................................

#..............................USE GSL..........................................
OPTION(USE_GSL "Activate GSL library" ON)

IF(USE_GSL) # set with: cmake . -DUSE_GSL:BOOL=ON
  MESSAGE("# Compiling with GSL Library...")
  FIND_PACKAGE (GSL)
  IF(GSL_FOUND)
    MESSAGE(STATUS "GSL found ----------")
    INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIR})
    INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
    LINK_DIRECTORIES(${GSL_LIBRARY_DIR})
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_GSL")
  ELSE()
    MESSAGE(STATUS "GSL NOT found ----------")
  ENDIF(GSL_FOUND)
ELSE()
  MESSAGE(STATUS "GSL is not activated ----------")
ENDIF(USE_GSL)
#...............................................................................

#..............................USE HF5..........................................
OPTION(USE_HDF5 "Compile with the HDF5 library" ON)
  
IF(USE_HDF5) # set with: cmake . -DUSE_HF5:BOOL=ON
  MESSAGE("Compiling with HDF5 Library")
  IF(HDF5_INCDIR)
    INCLUDE_DIRECTORIES(${HDF5_INCDIR})
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_HDF5") 
  ELSE()
    FIND_PACKAGE (HDF5)
    IF(HDF5_FOUND)
      MESSAGE(STATUS "HDF5 found -----------" ${HDF5_INCLUDE_DIRS})
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_HDF5")
      LIST(APPEND INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIRS})
      INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
    ELSE()
      MESSAGE(STATUS "HDF5 NOT found ----------")
    ENDIF(HDF5_FOUND)
  ENDIF(HDF5_INCDIR)
ENDIF(USE_HDF5)
#...............................................................................


#..............................USE MPI..........................................
OPTION(USE_MPI "Compile with the MPI library"  ON)

IF(USE_MPI)
  MESSAGE("Compiling with MPI Library")
  IF(MPI_INCDIR)
    INCLUDE_DIRECTORIES(${MPI_INCDIR})
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_MPI") 
  ELSE()
    FIND_PACKAGE (MPI REQUIRED)
    IF(MPI_FOUND)
      MESSAGE(STATUS "MPI found -----------" ${MPI_C_INCLUDE_PATH})
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_MPI")
      LIST(APPEND INCLUDE_DIRECTORIES ${MPI_C_INCLUDE_PATH})
      INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
    ELSE()
      MESSAGE(STATUS "MPI NOT found ----------")
    ENDIF(MPI_FOUND)
  ENDIF(MPI_INCDIR)
ENDIF(USE_MPI)
#...............................................................................                                        
# .......................................................
MESSAGE(STATUS "# Include directories found :" ${INCLUDE_DIRECTORIES})
# .......................................................

# .................ADDING OUR LIBRARIES..................
INCLUDE_DIRECTORIES(src)
# TODO: replace by the following line when necessary
# include_directories(src ../util/src)
# .......................................................


# ......................CREATING SCHNAPS LIB......................
SET(source_files "")

# We create a macro to add every *.c and *.h for a given library
MACRO(SCHNAPS_ADD_LIB _NAME)
  LIST(APPEND source_files src/${_NAME}.c)
  LIST(APPEND source_files src/${_NAME}.h)
ENDMACRO(SCHNAPS_ADD_LIB)

SCHNAPS_ADD_LIB(simulation)
SCHNAPS_ADD_LIB(simulation_spu)
SCHNAPS_ADD_LIB(linear_solver)
SCHNAPS_ADD_LIB(field)
SCHNAPS_ADD_LIB(interface)
SCHNAPS_ADD_LIB(implicit)
SCHNAPS_ADD_LIB(geometry)
SCHNAPS_ADD_LIB(global)
SCHNAPS_ADD_LIB(model)
SCHNAPS_ADD_LIB(interpolation)
SCHNAPS_ADD_LIB(macromesh)
SCHNAPS_ADD_LIB(maxwell)
SCHNAPS_ADD_LIB(waterwave2d)
SCHNAPS_ADD_LIB(mhd)
SCHNAPS_ADD_LIB(gyro)
SCHNAPS_ADD_LIB(pic)
SCHNAPS_ADD_LIB(skyline)
SCHNAPS_ADD_LIB(skyline_spu)
SCHNAPS_ADD_LIB(klu_csr)
#TODO: why are these out ?
#SCHNAPS_ADD_LIB(csparse)
#SCHNAPS_ADD_LIB(csparse_spu)
SCHNAPS_ADD_LIB(dpackfgmres)
SCHNAPS_ADD_LIB(advanced_linear_solver)
SCHNAPS_ADD_LIB(solvercontinuous)
SCHNAPS_ADD_LIB(solverpoisson)
SCHNAPS_ADD_LIB(quantities_vp)
SCHNAPS_ADD_LIB(diagnostics_vp)
SCHNAPS_ADD_LIB(collision)
SCHNAPS_ADD_LIB(lattice)
SCHNAPS_ADD_LIB(lbm_generic)
SCHNAPS_ADD_LIB(lbm_timeschemes)
SCHNAPS_ADD_LIB(lbm_diagnostics)
SCHNAPS_ADD_LIB(acoustic)
SCHNAPS_ADD_LIB(io)

LIST(APPEND
  source_files
  src/schnaps.h
  src/h20phi.h
  src/lbm_models_data.h
  )

# Adding SuiteSparse lib......
IF(SuiteSparse_FOUND)
  SCHNAPS_ADD_LIB(solverumfpack)
ENDIF(SuiteSparse_FOUND)

# Adding OpenCl libs......
IF (OPENCL_FOUND)
  SCHNAPS_ADD_LIB(field_cl)
  SCHNAPS_ADD_LIB(clutils)
  SCHNAPS_ADD_LIB(clinfo)
ENDIF (OPENCL_FOUND)

IF(PYTHONLIBS_FOUND)
    SCHNAPS_ADD_LIB(python_call)
ENDIF(PYTHONLIBS_FOUND)

#................Creating schnaps library.................
ADD_LIBRARY (_schnaps ${source_files})
# ........................................................

#................Adding external librares.................
SET(LIBS "")



IF (PYTHONLIBS_FOUND)
  LIST(APPEND LIBS _schnaps ${PYTHON_LIBRARY} m)
ELSE()
  LIST(APPEND LIBS _schnaps  m)
ENDIF (PYTHONLIBS_FOUND)

IF (OPENCL_FOUND)
  LIST(APPEND LIBS _schnaps ${OPENCL_LIBRARIES} m)
ELSE()
  LIST(APPEND LIBS _schnaps  m)
ENDIF (OPENCL_FOUND)

IF (STARPU_FOUND)
  LIST(APPEND LIBS ${STARPU_SHM_LIBRARIES})
  LIST(APPEND LIBS ${STARPU_LIBRARIES})
  MESSAGE(STATUS "# StarPU SHM found at:" ${STARPU_SHM_LIBRARIES})
ENDIF (STARPU_FOUND)

IF (HAVE_IGRAPH)
  LIST(APPEND LIBS  ${IGRAPH_LIBRARIES})
ENDIF (HAVE_IGRAPH)
IF (GSL_FOUND)
  LIST(APPEND LIBS ${GSL_LIBRARIES})
ENDIF (GSL_FOUND)

IF (HDF5_FOUND)
    MESSAGE(" HDF5 libs found at :" ${HDF5_LIBRARIES})
  LIST(APPEND LIBS ${HDF5_LIBRARIES})
ENDIF (HDF5_FOUND)
# nearest neighbour fast search library
FIND_PACKAGE(flann)
IF(FLANN_FOUND)
  MESSAGE(STATUS "# Delicious!")
  LIST(APPEND LIBS ${FLANN_LIBRARY})
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_FLANN")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -D_WITH_FLANN")
ENDIF()

IF(SuiteSparse_FOUND)
  LIST(APPEND LIBS ${SuiteSparse_LIBRARIES})
  LIST(APPEND LIBS ${SuiteSparse_LAPACK_BLAS_LIBRARIES})
ENDIF()

# LIST(APPEND LIBS "/home/helluy/schnaps/old_schnaps/build/model_trans.so")
# LINK_DIRECTORIES("/home/helluy/schnaps/old_schnaps/build/")
# .......................................................

# .......................................................
ENABLE_TESTING()
# .......................................................

# .....................ADDING MAIN EXEC..................
# main program schnaps
ADD_EXECUTABLE(schnaps src/schnaps.c)
TARGET_LINK_LIBRARIES(schnaps ${LIBS})
# .......................................................

# .................ADDING OTHER TESTS....................
# Find all of the files ending in .c in the test directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "test/")
file(GLOB globtest test/*.c)
MESSAGE(STATUS "globtest")
foreach(test IN LISTS globtest)
  string(REPLACE  ".c"  ""  test  ${test})
  string(REGEX REPLACE  "([^/]+$)"  ""  dir  ${test})
  string(REPLACE  ${dir}  ""  test  ${test})
  #MESSAGE(STATUS ${test})
  list(APPEND testlist ${test})
endforeach(test)

foreach(test IN LISTS testlist)
  MESSAGE(STATUS "Adding test " ${test})
  add_executable(
    ${test}
    test/test.h
    test/${test}.c
    )
  target_link_libraries(${test} ${LIBS})
  add_test(test/${test} test/${test})
  set_property(TEST test/${test} PROPERTY RUN_SERIAL 1)
endforeach(test)

# .......................ADDING EXAMPLES....................
# Find all of the files ending in .c in the example directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "example/")
file(GLOB globexample example/*.c)
MESSAGE(STATUS "globexample")
foreach(example IN LISTS globexample)
  string(REPLACE  ".c"  ""  example  ${example})
  string(REGEX REPLACE  "([^/]+$)"  ""  dir  ${example})
  string(REPLACE  ${dir}  ""  example  ${example})
  #MESSAGE(STATUS ${example})
  list(APPEND examplelist ${example})
endforeach(example)

foreach(example IN LISTS examplelist)
  MESSAGE(STATUS "Adding example " ${example})
  add_executable(
    ${example}
    test/test.h
    example/${example}.c
    )
  target_link_libraries(${example} ${LIBS})
endforeach(example)
# .......................................................

#...............ADDING RELEASE/DEBUG FLAGS...............
SET(CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS}  -Wall")
SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} ")
# .......................................................

#................INSTALLING OPENCL FILES.................
FILE(GLOB INPUT_FILES *.sh)
FILE(COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
#.......................................................


MESSAGE("")
MESSAGE(STATUS "############################################################")
MESSAGE(STATUS "####")
MESSAGE(STATUS "#### CMAKE_BUILD_TYPE    : ${CMAKE_BUILD_TYPE}")
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
   MESSAGE(STATUS "#### ${CMAKE_C_FLAGS_DEBUG}")
ELSE()
   MESSAGE(STATUS "#### ${CMAKE_C_FLAGS}")
ENDIF()
MESSAGE(STATUS "#### PythonLibs Found    : ${PYTHONLIBS_FOUND}")
MESSAGE(STATUS "#### SuiteSparse Found   : ${SuiteSparse_FOUND}")
MESSAGE(STATUS "#### OpenCl Found        : ${OPENCL_FOUND}")
MESSAGE(STATUS "#### StarPU Found        : ${STARPU_FOUND}")
MESSAGE(STATUS "#### igraph Found        : ${HAVE_IGRAPH}")
MESSAGE(STATUS "#### GC USE              : ${GC_ENABLED}")
IF(USE_GSL)
MESSAGE(STATUS "#### GSL Found           : ${GSL_FOUND}")
ELSE()
MESSAGE(STATUS "#### GSL USE             : ${USE_GSL}")
ENDIF()
MESSAGE(STATUS "#### HDF5 Found          : ${HDF5_FOUND}") 
MESSAGE(STATUS "#### MPI Found           : ${MPI_FOUND}")                                           
MESSAGE(STATUS "####")
MESSAGE(STATUS "############################################################")
