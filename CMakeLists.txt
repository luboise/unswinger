cmake_minimum_required(VERSION 3.28)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_CXX_EXTENSIONS Off)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Build Configurations" FORCE)

project(MusicUnswingerProject)

add_executable(MusicUnswinger)

# libsndfile
find_package(SndFile REQUIRED)
include_directories(${SNDFILE_INCLUDE_DIRS})

if(WIN32)
	target_link_directories(MusicUnswinger PRIVATE ext/FFTW)
	target_include_directories(MusicUnswinger PRIVATE ext/FFTW)
	target_link_libraries(MusicUnswinger PRIVATE libfftw3-3.lib)
elseif(UNIX)
	target_link_libraries(MusicUnswinger PRIVATE fftw3)
endif()

target_link_libraries(MusicUnswinger
	PUBLIC	${SNDFILE_LIBRARIES}
)

add_subdirectory(src)