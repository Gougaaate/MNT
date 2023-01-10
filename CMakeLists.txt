cmake_minimum_required(VERSION 3.10)
project(create_raster)

# Add source files
set(SOURCES
        src/main.cpp
        src/Point.cpp
        src/projection.cpp
        src/hillshading.cpp
        src/display.cpp
        )

# Add header files
set(HEADERS
        src/Point.hpp
        src/delaunator.hpp
        src/hillshading.hpp
        src/projection.hpp
        src/display.hpp
        )

# Add executable
add_executable(create_raster ${SOURCES} ${HEADERS})

target_link_libraries(create_raster proj)