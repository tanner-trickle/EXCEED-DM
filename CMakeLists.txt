cmake_minimum_required(VERSION 3.16)

# Project information ###########################

project(EXCEED-DM
        VERSION 0.2.6
        DESCRIPTION "OpenMPI Fortran program that computes Dark Matter (DM) induced electronic transition rates for any gapped crystal target."
        HOMEPAGE_URL "https://tanner-trickle.github.io/EXCEED-DM/index.html"
        LANGUAGES Fortran)

#################################################

# Build type options ############################

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -Wall -Wextra")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")

#################################################

# Main executable ###############################

file(GLOB_RECURSE sources src/*.f90)
add_executable(exdm ${sources})

# External Libraries ############################

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMakeModules/")

# Automatically find required libraries.
find_package(FFTW REQUIRED)
find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS Fortran HL)
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

set(ext_libraries 
    ${FFTW_LIBRARIES}
    ${MPI_Fortran_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
)

set(ext_include_dirs
    ${FFTW_INCLUDE_DIRS}
    ${MPI_Fortran_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
)

# Link all the libraries.
target_link_libraries(exdm PUBLIC ${ext_libraries})

# Location of .mod files
target_include_directories(exdm PUBLIC ${ext_include_dirs})

#################################################
