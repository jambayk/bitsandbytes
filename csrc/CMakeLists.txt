cmake_minimum_required(VERSION 3.12)

# Include the ExternalProject module
include(ExternalProject)

if(CMAKE_BUILD_TYPE MATCHES Debug)
    set(OUTPUT_SUFFIX "d")
endif()

if(WIN32)
# Define a target for the pthreads library
ExternalProject_Add(
    pthreads
    PREFIX ${CMAKE_BINARY_DIR}/externals/pthreads
    GIT_REPOSITORY https://github.com/GerHobbelt/pthread-win32.git
    GIT_TAG version-3.1.0-release  # or specify a specific release/tag
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -DCMAKE_BUILD_TYPE=Release
)
endif()

# Define a target for the CUB library
ExternalProject_Add(
    cub_library
    PREFIX ${CMAKE_BINARY_DIR}/externals/cub
    GIT_REPOSITORY https://github.com/NVIDIA/cub.git
    GIT_TAG 2.0.1  # specify the desired release/tag
    CONFIGURE_COMMAND ""  # Skip configure step
    BUILD_COMMAND ""      # Skip build step
    INSTALL_COMMAND ""    # Skip install step
    LOG_DOWNLOAD ON       # Log the download step for debugging
)

add_library(bitsandbytes SHARED
    ops.cu
    kernels.cu
    common.cpp
    cpu_ops.cpp
    pythonInterface.cpp
)

# Ensure CUDA toolkit is available
set(CMAKE_VS_PLATFORM_TOOLSET_CUDA "12" CACHE STRING "CUDA platform toolset version")
if(WIN32)	
else()	
	set(CUDA_TOOLKIT_ROOT_DIR $ENV{CUDA_HOME})
endif()

# Set a flag to indicate whether cublaslt is included
set(USE_CUBLASLT TRUE CACHE BOOL "Enable/Disable cublasLT")  # Change to OFF if you want to disable cublaslt

# Conditional inclusion of cublaslt
if(USE_CUBLASLT)
    target_compile_definitions(bitsandbytes PRIVATE USE_CUBLASLT)
    target_link_libraries(bitsandbytes PRIVATE cublasLt)
else()
    target_compile_definitions(bitsandbytes PRIVATE NO_CUBLASLT)
endif()

# Remove dots from the CUDA version
string(REPLACE "." "" CMAKE_VS_PLATFORM_TOOLSET_CUDA_NO_DOTS ${CMAKE_VS_PLATFORM_TOOLSET_CUDA})

# Set the output name for the shared library without dots
set(OUTPUT_NAME_SUFFIX "_nocublaslt")
if(USE_CUBLASLT)
    set(OUTPUT_NAME_SUFFIX "")
endif()

set_target_properties(bitsandbytes PROPERTIES
    OUTPUT_NAME "bitsandbytes_cuda${CMAKE_VS_PLATFORM_TOOLSET_CUDA_NO_DOTS}${OUTPUT_NAME_SUFFIX}"
)

# Include directories
target_include_directories(bitsandbytes PUBLIC
    ${CUDA_TOOLKIT_ROOT_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_BINARY_DIR}/externals/cub/src/cub_library  # Include CUB headers
    ${CMAKE_BINARY_DIR}/externals/pthreads/include  # Include pthreads
)

# Link directories
if(WIN32)
    target_link_directories(bitsandbytes PUBLIC 
        ${CUDA_TOOLKIT_ROOT_DIR}/lib64
        ${CMAKE_BINARY_DIR}/externals/pthreads/lib  # Link pthreads
    )	
else()
    target_link_directories(bitsandbytes PUBLIC 
        ${CUDA_TOOLKIT_ROOT_DIR}/lib64
    )
endif()

# Compiler definitions
target_compile_definitions(bitsandbytes PRIVATE BUILD_CUDA)

# Link libraries
target_link_libraries(bitsandbytes PRIVATE
    ${CUDA_LIBRARIES}
    cudart cublas curand cusparse
)

# Link with pthreads
if(WIN32)
	if(CMAKE_BUILD_TYPE MATCHES Debug)
		target_link_libraries(bitsandbytes PRIVATE pthreadVC3d)
	else()
		target_link_libraries(bitsandbytes PRIVATE pthreadVC3)
	endif()
else()
    target_link_libraries(bitsandbytes PRIVATE pthread)
endif()

if(MSVC)
    set_property(TARGET bitsandbytes APPEND PROPERTY LINK_OPTIONS /OPT:REF)
endif()

set_target_properties(bitsandbytes
    PROPERTIES 
    CUDA_SEPARABLE_COMPILATION ON 
    CUDA_ARCHITECTURES all
)
