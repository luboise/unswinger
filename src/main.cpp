#include <iostream>

// ./unswinger filepath bpm offset
#define NUM_ARGS 3

void printHelp(void) {
    std::cerr << "./unswinger filepath bpm offset" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != NUM_ARGS) {
        printHelp();
        return EXIT_FAILURE;
    }

        for (size_t i = 0; i < argc; i++) {
        auto arg = argv[i];
    }

    return EXIT_SUCCESS;
}