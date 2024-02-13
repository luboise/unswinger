#include <iostream>
#include <string>

#include <filesystem>
namespace fs = std::filesystem;

#include <sound.h>

#include <fftw3.h>

// ./unswinger filepath bpm offset
#define NUM_USER_ARGS 4
#define NUM_FFT_ARGS 3

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
    if (argc != NUM_USER_ARGS + 1 && argc != NUM_USER_ARGS + NUM_FFT_ARGS + 1) {
        printHelp();
        return EXIT_FAILURE;
    }

    std::string inpath;

    std::string removeStr;
    bool removeSwing;

    double songBPM;
    double offset;

    bool userSpecifiedFft = (argc == NUM_USER_ARGS + NUM_FFT_ARGS + 1);
    SoundFile::FFTParams fftParams;

    try {
        inpath = argv[1];

        removeStr = std::string(argv[2]);
        removeSwing = removeStr == "remove";

        songBPM = std::stod(argv[3]);
        offset = std::stod(argv[4]);

        if (argc == NUM_USER_ARGS + NUM_FFT_ARGS + 1) {
            fftParams.fftWindowSize = std::stoi(argv[NUM_USER_ARGS + 1]);
            fftParams.overlayWhenFast = std::stod(argv[NUM_USER_ARGS + 2]);
            fftParams.overlayWhenSlow = std::stod(argv[NUM_USER_ARGS + 3]);
        }

    } catch (std::exception& e) {
        std::cerr << "Failed to process user arguments with error: \n"
                  << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Loading file: " << inpath << std::endl;
    SoundFile file(inpath);
    if (userSpecifiedFft) file.setFftParams(fftParams);

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

    file.addSwingVocoded(songBPM, offset, removeSwing);

    file.normalise();

    auto inpathFS = fs::path(inpath);

    std::string fileStem = inpathFS.stem().string() + "_swung";
    if (userSpecifiedFft) {
        std::stringstream ss;
        ss << "_" << fftParams.fftWindowSize << "_" << fftParams.overlayWhenFast
           << "_"
           << fftParams.overlayWhenSlow;
        fileStem += ss.str();
    }

    auto outpath = fs::path(fileStem + inpathFS.extension().string());

    file.exportToFile(outpath.string());

    return EXIT_SUCCESS;
}