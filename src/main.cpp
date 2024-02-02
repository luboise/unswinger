#include <iostream>
#include <string>

#include <sound.h>

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

    std::string filepath;
    double songBPM;
    double offset;

    try {
        filepath = argv[1];
        songBPM = std::stod(argv[2]);
        offset = std::stod(argv[3]);
    } catch (std::exception& e) {
        std::cerr << "Failed to process user arguments with error: \n"
                  << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    SoundFile file(filepath);
    file.addSwing(songBPM, offset);
    file.exportToFile("out.flac");

    return EXIT_SUCCESS;
}