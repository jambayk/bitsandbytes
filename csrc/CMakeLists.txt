cmake_minimum_required(VERSION 3.12)

# Ensure CUDA toolkit is available
set(CMAKE_VS_PLATFORM_TOOLSET_CUDA "12" CACHE STRING "CUDA platform toolset version")

# Include the ExternalProject module
include(ExternalProject)

# Define a target for the pthreads library
ExternalProject_Add(
    pthreads
    PREFIX ${CMAKE_BINARY_DIR}/externals/pthreads
    GIT_REPOSITORY https://github.com/GerHobbelt/pthread-win32.git
    GIT_TAG master  # or specify a specific release/tag
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -DCMAKE_FIND_PACKAGE_PreferConfig=TRUE
)

add_library(bitsandbytes SHARED
    ops.cu
    kernels.cu
    common.cpp
    cpu_ops.cpp
    pythonInterface.cpp
)

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
    cudart cublas cublasLt curand cusparse
)

# Link with pthreads
target_link_libraries(bitsandbytes PRIVATE pthreadVC3)

if(MSVC)
    set_property(TARGET bitsandbytes APPEND PROPERTY LINK_OPTIONS /OPT:REF)
endif()

set_target_properties(bitsandbytes
    PROPERTIES 
    CUDA_SEPARABLE_COMPILATION ON 
    CUDA_ARCHITECTURES all
)
