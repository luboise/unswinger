cmake_minimum_required(VERSION 3.28)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_CXX_EXTENSIONS Off)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Build Configurations" FORCE)

project(MusicUnswingerProject)

add_executable(MusicUnswinger)

# Libsndfile
add_compile_definitions(ENABLE_EXTERNAL_LIBS ENABLE_MPEG)
add_subdirectory(ext/libsndfile)

target_link_libraries(MusicUnswinger
	PUBLIC	libsndfile
)

add_subdirectory(src)