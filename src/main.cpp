#include <iostream>
#include <string>

#include <filesystem>
namespace fs = std::filesystem;

#include <sound.h>

#include <fftw3.h>

// ./unswinger filepath bpm offset
#define NUM_USER_ARGS 4

void printHelp(void) {
    std::cerr << "./unswinger filepath status bpm offset"
                 "\nfilepath - Path to the audio you want to add swing to"
                 "\nstatus - Add or remove, add if adding swing, remove if "
                 "removing swing"
                 "\nbpm - The Beats Per Minute of the audio"
                 "\noffset - The sync offset of the bpm (in seconds)"
              << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != NUM_USER_ARGS + 1) {
        printHelp();
        return EXIT_FAILURE;
    }

    std::string inpath;

    std::string removeStr;
    bool removeSwing;

    double songBPM;
    double offset;

    try {
        inpath = argv[1];

        removeStr = std::string(argv[2]);
        removeSwing = removeStr == "remove";

        songBPM = std::stod(argv[3]);
        offset = std::stod(argv[4]);
    } catch (std::exception& e) {
        std::cerr << "Failed to process user arguments with error: \n"
                  << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    SoundFile file(inpath);

    file.addSwingFourier(songBPM, offset, removeSwing);
    file.normalise();
    // file.setChannel(0, 0, file.getPitched(file.getChannel(0), 5));
    // file.setChannel(1, 0, file.getPitched(file.getChannel(1), 5));

    auto inpathFS = fs::path(inpath);

    auto outpath = fs::path(inpathFS.stem().string() + "_swung" +
                            inpathFS.extension().string());

    file.exportToFile(outpath.string());

    return EXIT_SUCCESS;
}