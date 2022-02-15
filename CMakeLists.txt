cmake_minimum_required(VERSION 3.16)

# Project information ###########################

project(EXCEED-DM
        VERSION 0.3.0
        DESCRIPTION "OpenMPI Fortran program to compute Dark Matter (DM) induced electronic transition rates in gapped crystal targets."
        HOMEPAGE_URL "https://exceed-dm.caltech.edu"
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

set(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_SOURCE_DIR}/build/mod")

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
