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

    std::cout << "Loading file: " << inpath << std::endl;
    SoundFile file(inpath);

    if (!file.isValid()) {
        std::cerr << "Unable to proceed." << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Successful file load.\n"
              << "Audio length: " << file.getSoundLength()
              << "\nTotal sample count: " << file.getSampleCount()
              << "\nChannels: " << file.getChannelCount() << " ("
              << file.getSampleCount() / file.getChannelCount()
              << " samples per channel)" << std::endl;

    //file.addSwingVocoded(songBPM, offset, removeSwing);

    for (size_t channel = 0; channel < file.getChannelCount(); channel++)
        file.setChannel(channel, 0, file.getVocoded(file.getChannel(channel), 1.5));

    file.normalise();

    auto inpathFS = fs::path(inpath);

    auto outpath = fs::path(inpathFS.stem().string() + "_swung" +
                            inpathFS.extension().string());

    file.exportToFile(outpath.string());

    return EXIT_SUCCESS;
}