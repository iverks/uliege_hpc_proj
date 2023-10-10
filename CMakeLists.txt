cmake_minimum_required(VERSION 3.5)

# Definitions
set(PROJECT_NAME fdtd)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(${PROJECT_NAME})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O3 -Wall -Wextra")

# Include paths
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/)

# Assign the include directories
include_directories(${INCLUDE_DIRS})

# Create the executable
set(SOURCE_FILES fdtd.c)
add_executable(fdtd ${SOURCE_FILES})
target_link_libraries(fdtd m)
