#include <iostream>
#include <string>

#include <filesystem>
namespace fs = std::filesystem;

#include <sound.h>

#include <fftw3.h>

// ./unswinger filepath bpm offset
#define NUM_USER_ARGS 3

void printHelp(void) {
    std::cerr << "./unswinger filepath bpm offset\n"
                 "filepath - Path to the audio you want to add swing to\n"
                 "bpm - The Beats Per Minute of the audio\n"
                 "offset - The sync offset of the bpm (in seconds)"
              << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != NUM_USER_ARGS + 1) {
        printHelp();
        return EXIT_FAILURE;
    }

    std::string inpath;
    double songBPM;
    double offset;

    try {
        inpath = argv[1];
        songBPM = std::stod(argv[2]);
        offset = std::stod(argv[3]);
    } catch (std::exception& e) {
        std::cerr << "Failed to process user arguments with error: \n"
                  << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    SoundFile file(inpath);

    file.addSwingFourier(songBPM, offset);
    file.normalise();
    //file.setChannel(0, 0, file.getPitched(file.getChannel(0), 5));
    //file.setChannel(1, 0, file.getPitched(file.getChannel(1), 5));

    auto outpath =
        fs::path("out").replace_extension(fs::path(inpath).extension());

    file.exportToFile(outpath.string());

    return EXIT_SUCCESS;
}