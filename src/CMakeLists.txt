target_sources(MusicUnswinger PRIVATE main.cpp sound.h sound.cpp utils.h utils.cpp fft/fftplan.cpp fft/fftplan.h)

target_include_directories(MusicUnswinger
PRIVATE	${CMAKE_CURRENT_SOURCE_DIR}
)