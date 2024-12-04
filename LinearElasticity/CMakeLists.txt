## CMAKE for Dendro-KT
cmake_minimum_required(VERSION 2.8)
project(le-kt)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/dendrite-kt/cmake-modules")
# For now we just make it compulsory to have LAPACK installed.
#Later we will make it possible if LAPACK is not present to automaticall install before compiling dendro5
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()


if (MPI_COMPILE_FLAGS)
    set(COMPILE_FLAGS "${COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
endif ()

if (MPI_LINK_FLAGS)
    set(LINK_FLAGS "${LINK_FLAGS} ${MPI_LINK_FLAGS}")
endif ()


## options for dendro

option(USE_64BIT_INDICES "Use 64-Bit indices. Reverts to 32-bit if turned off" ON)
option(ALLTOALLV_FIX "Use K-way all to all v" ON)
option(SPLITTER_SELECTION_FIX "Turn on Splitter Selection fix" ON)
option(DIM_2 "use the two dimentional sorting" OFF)
option(WITH_BLAS_LAPACK "build using BLAS and LAPACk" ON)
option(MANUAL_BLAS_LAPACK "configure BLAS and LAPACK Manually" OFF)
option(DENDRO_VTK_BINARY "write vtk/vtu files in binary mode " ON)
option(DENDRITE_VTU_ASCII "write vtk/vtu files in ASCII mode " OFF)
option(DENDRO_VTK_ZLIB_COMPRES "write vtk/vtu files in binary mode with zlib compression (only compatible with binary mode) " OFF)
option(BUILD_WITH_PETSC " build dendro with PETSC " ON)
option(HILBERT_ORDERING "use the Hilbert space-filling curve to order orthants" OFF)
option(BUILD_EXAMPLES "build example programs" ON)
option(ENABLE_4D "enable 4D computation" OFF)
option(ENABLE_2D "enable 2D computation" OFF)
option(ENABLE_3D "enable 3D computation" ON)
option(TENSOR "Use Tensor Operation" OFF)
option(PROFILING "Enable profiling" OFF)
option(IBM "Enable IBM Functions" OFF)

set(KWAY 128 CACHE INT 128)
set(NUM_NPES_THRESHOLD 2 CACHE INT 2)

#set the build type to release by default.
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
            "Choose the type of build, options are: Debug Release " FORCE)
endif ()

if (WITH_BLAS_LAPACK)
    add_definitions(-DWITH_BLAS_LAPACK)
    if (DEFINED ENV{MKLROOT})
        find_package(LAPACK COMPONENTS MKL REQUIRED)
        set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
        if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mkl")
        else ()
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm")
        endif ()
        message(STATUS ${LAPACK_LIBRARIES})
    elseif (MANUAL_BLAS_LAPACK)
        if ("$ENV{BLAS}" STREQUAL "")
            message("Environment Variable BLAS is not set. Please set it to BLAS directory")
        endif ()
        if ("$ENV{LAPACK}" STREQUAL "")
            message("Enviroment Variable LAPACK is note set. Please set it to LAPACK directory. ")
        endif ()
        set(LAPACKE_DIR $ENV{LAPACK}/LAPACKE)
        set(BLAS_LIBS $ENV{BLAS}/lib)
        set(LAPACK_LIBS $ENV{LAPACK}/lib)
        set(LAPACK_LINKER_FLAGS -llapacke -llapack -lblas -lgfortran -lquadmath)
        set(LAPACK_LIBRARIES ${LAPACK_LIBS}/liblapacke.a ${LAPACK_LIBS}/liblapack.a ${BLAS_LIBS}/libblas.a -static libgfortran.a libquadmath.a)
        set(LINK_FLAGS "${LINK_FLAGS} ${LAPACK_LINKER_FLAGS}")
    else ()
        find_package(BLAS REQUIRED)
        find_package(LAPACK REQUIRED)
        set(LAPACK_LINKER_FLAGS -llapacke -llapack -lblas -lgfortran -lquadmath)
        set(LAPACKE_DIR $ENV{LAPACK}/LAPACKE)
        set(LINK_FLAGS "${LINK_FLAGS} ${LAPACK_LINKER_FLAGS}")
        find_library(LAPACKE_LIB
                NAMES lapacke lapackelib liblapacke
                HINTS "/usr/lib/"
                )
        set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${LAPACKE_LIB})
        message(STATUS ${LAPACK_LIBRARIES})
    endif ()
endif ()

if (BUILD_WITH_PETSC)
    find_package(PETSc REQUIRED)
    add_definitions(-DBUILD_WITH_PETSC)
endif ()

if (TENSOR)
    add_definitions(-DTENSOR)
    message("Enabling Tensor operation")
endif ()
if (DIM_2)
    add_definitions(-DDIM_2)
endif ()

if (PROFILING)
    add_definitions(-DPROFILING)
endif ()
if (USE_64BIT_INDICES)
    add_definitions(-DUSE_64BIT_INDICES)
    #message('Configuring 64BIT indices')
endif ()

if (ALLTOALLV_FIX)
    add_definitions(-DALLTOALLV_FIX)
    add_definitions(-DKWAY=${KWAY})
endif ()

if (SPLITTER_SELECTION_FIX)
    add_definitions(-DSPLITTER_SELECTION_FIX)
    add_definitions(-DNUM_NPES_THRESHOLD=${NUM_NPES_THRESHOLD})
endif ()

if (ALLTOALL_SPARSE)
    add_definitions(-DALLTOALL_SPARSE)
endif ()

if (DENDRO_VTK_BINARY)
else ()
    set(DENDRO_VTK_ZLIB_COMPRES OFF)
endif ()

if (DENDRITE_VTU_ASCII)
    add_definitions(-DDENDRITE_VTU_ASCII)
    message("Writing in ASCII format")
else ()
    message("Writing in Binary format")
endif ()


if (DENDRO_VTK_BINARY)
    add_definitions(-DDENDRO_VTU_BINARY)
    if (DENDRO_VTK_ZLIB_COMPRES)
        add_definitions(-DDENDRO_VTU_ZLIB)
    endif ()
else ()
    add_definitions(-DDENDRO_VTU_ASCII)
endif ()

if (ENABLE_4D)
    add_definitions(-DENABLE_4D)
    message("Enabling 4D computation")
    set(ENABLE_4D,ON)
    set(ENABLE_3D,OFF)
elseif (ENABLE_2D)
    add_definitions(-DENABLE_2D)
    message("Enabling 2D computation")
    set(ENABLE_2D,ON)
    set(ENABLE_3D,OFF)
elseif (ENABLE_3D)
    add_definitions(-DENABLE_3D)
    set(ENABLE_3D,ON)
    message("Enabling 3D computation")
endif ()

if (HILBERT_ORDERING)
    add_definitions(-DHILBERT_ORDERING)
endif ()

if (IBM)
    add_definitions(-DIBM)
    message("Enabling IMGA Computation")
endif ()

set(DENDRITEkT_BUILD_EXAMPLES OFF CACHE BOOL "Build Dendrite examples")
add_subdirectory(dendrite-kt "${CMAKE_CURRENT_BINARY_DIR}/dendrite-kt")

set(LE_INC
        include/util.h
        include/LEInputData.h
        include/LENodeData.h
        include/LEBCSetup.h
        include/TSLEEquation.h
        include/SSLEEquation.h
        include/LERefine.h
        include/LEInputDataStructs.h
        )

set(TLE_SRC src/TSLE.cpp)

add_executable(tle-kt ${TLE_SRC} ${LE_INC})
target_include_directories(tle-kt PUBLIC include)
target_link_libraries(tle-kt dendriteKT dendroKT ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)


set(SLE_SRC src/SSLE.cpp include/SBMCalc.h include/GetSurfaceGp.h include/GetTrueSurfaceGP.h include/BFS.h include/ElementMarker.h include/CheckSurface.h
        include/CalcStress.h
        include/DARefine.h
        include/DACoarse.h)

add_executable(sle-kt ${SLE_SRC} ${LE_INC})
target_include_directories(sle-kt PUBLIC include)
target_link_libraries(sle-kt dendriteKT dendroKT ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
