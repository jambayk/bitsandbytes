cmake_minimum_required(VERSION 3.12)

if(NOT WIN32)
    message(FATAL_ERROR "This CMake project is intended for Windows only! for Linux use the makefile.")
endif()

# Ensure CUDA toolkit is available
set(CMAKE_VS_PLATFORM_TOOLSET_CUDA "12" CACHE STRING "CUDA platform toolset version")

# Set a flag to indicate whether cublaslt is included
set(USE_CUBLASLT OFF CACHE BOOL "Enable/Disable cublasLT")  # Change to OFF if you want to disable cublaslt

# Include the ExternalProject module
include(ExternalProject)

# Define a target for the pthreads library
ExternalProject_Add(
    pthreads
    PREFIX ${CMAKE_BINARY_DIR}/externals/pthreads
    GIT_REPOSITORY https://github.com/GerHobbelt/pthread-win32.git
    GIT_TAG a301ede8e6a285fe0891c81eaba847d1b702c313  # or specify a specific release/tag
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
)

add_library(bitsandbytes SHARED
    ops.cu
    kernels.cu
    common.cpp
    cpu_ops.cpp
    pythonInterface.c
)

SET_SOURCE_FILES_PROPERTIES(pythonInterface.c PROPERTIES LANGUAGE CXX)

# Include directories
target_include_directories(bitsandbytes PUBLIC
    ${CUDA_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/cub
    ${CMAKE_BINARY_DIR}/externals/pthreads/include  # Include pthreads
)

# Link directories
target_link_directories(bitsandbytes PUBLIC 
	${CUDA_TOOLKIT_ROOT_DIR}/lib64
	${CMAKE_BINARY_DIR}/externals/pthreads/lib  # Link pthreads
)

# Compiler definitions
target_compile_definitions(bitsandbytes PRIVATE BUILD_CUDA)

# Link libraries
target_link_libraries(bitsandbytes PRIVATE
    ${CUDA_LIBRARIES}
    cudart cublas curand cusparse
)

# Conditional inclusion of cublaslt
if(USE_CUBLASLT)
    target_compile_definitions(bitsandbytes PRIVATE USE_CUBLASLT)
    target_link_libraries(bitsandbytes PRIVATE cublasLt)
else()
    target_compile_definitions(bitsandbytes PRIVATE NO_CUBLASLT)
endif()


# Set the output name for the shared library without dots
set(OUTPUT_NAME_SUFFIX "_nocublaslt")
if(USE_CUBLASLT)
    set(OUTPUT_NAME_SUFFIX "")
endif()

# Link with pthreads
target_link_libraries(bitsandbytes PRIVATE pthreadVC3)

# Remove dots from the CUDA version
string(REPLACE "." "" CMAKE_VS_PLATFORM_TOOLSET_CUDA_NO_DOTS ${CMAKE_VS_PLATFORM_TOOLSET_CUDA})

# TODO: maybe allow the user to specify cuda architectures
# all does not work on old machine with GTX 1660 Ti 
# had to use 75
set_target_properties(bitsandbytes
    PROPERTIES 
    RUNTIME_OUTPUT_DIRECTORY "$<1:${CMAKE_SOURCE_DIR}/bitsandbytes>"
    OUTPUT_NAME "libbitsandbytes_cuda${CMAKE_VS_PLATFORM_TOOLSET_CUDA_NO_DOTS}${OUTPUT_NAME_SUFFIX}"
    CUDA_SEPARABLE_COMPILATION ON 
    CUDA_ARCHITECTURES all
)