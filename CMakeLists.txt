cmake_minimum_required(VERSION 3.7)
project(2 Fortran)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(OpenMP)
if(OpenMP_Fortran_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES INTEL)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl")
endif(CMAKE_Fortran_COMPILER_ID MATCHES INTEL)


if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fcheck=all")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -march=native")
    message(WARNING "Assuming cmake is run with FC=/path/to/caf")
endif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/regression)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/neuralnetwork)

add_library(ising_basis ${CMAKE_CURRENT_SOURCE_DIR}/lib/ising_basis.f90)
add_library(ising ${CMAKE_CURRENT_SOURCE_DIR}/lib/ising.f90)

target_link_libraries(ising_basis regression)

link_directories(/usit/abel/u1/anjohan/OpenCoarrays/build/lib64)
set(libs ising caf_mpi)

file(GLOB program_files ${CMAKE_CURRENT_SOURCE_DIR}/programs/*.*90)

foreach(program_file ${program_files})
    get_filename_component(program ${program_file} NAME_WE)
    add_executable(${program} ${program_file})
    target_link_libraries(${program} "${libs}")
endforeach(program_file)

target_link_libraries(linreg ising_basis)
target_link_libraries(logreg regression)
foreach(prog reg_nn_test_couplings reg_nn_convergence reg_nn_test_spins)
    target_link_libraries(${prog} neuralnetwork)
endforeach(prog)
