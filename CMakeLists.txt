cmake_minimum_required(VERSION 3.10)
project(create_raster)

# Add source files
set(SOURCES
        src/main.cpp
        src/Point.cpp
        )

# Add header files
set(HEADERS
        src/Point.h
        src/delaunator.hpp
        )

# Add executable
add_executable(create_raster ${SOURCES} ${HEADERS})

target_link_libraries(create_raster proj)