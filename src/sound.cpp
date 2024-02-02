#include "sound.h"

#include <iostream>

#include <filesystem>

namespace fs = std::filesystem;

SoundFile::SoundFile(const std::string& filepath) {
    if (!fs::exists(filepath) || fs::is_directory(filepath)) {
        std::cerr << "Unable to locate file: \"" << filepath
                  << "\" or it is a directory and not a file." << std::endl;
        return;
    }

    SF_INFO info;

    SNDFILE* snd = sf_open(filepath.c_str(), SFM_READ, &info);
    if (snd) {
        const uint32_t sample_count = info.frames * info.channels;
        std::vector<double> samples(sample_count);

        auto read_count = sf_readf_double(snd, samples.data(), info.frames);
        if (read_count != info.frames) {
            std::cerr << "Read failed for file \"" << filepath
                      << "\", frames read: " << read_count << "/" << info.frames
                      << std::endl;
            return;
        }

        this->_sndinfo = info;
        this->_samples = std::move(samples);

        sf_close(snd);

    } else {
        std::cerr << "Error opening sound file \"" << filepath << "\""
                  << std::endl;
    }
}

SoundFile::~SoundFile() {}

inline bool SoundFile::isValid() const { return this->_samples.size() > 0; }

void SoundFile::addSwing(const double bpm, double offset) {
    double beat_length;

    beat_length = 60 / bpm;

    if (offset > 0) {
        offset = std::fmod(offset, beat_length) - beat_length;
    }

    uint32_t beat = 0;

    // Trackers in seconds
    double tracker_l = offset;
    double tracker_r = -1;

    uint32_t left_frame = -1;
    uint32_t right_frame = -1;

    const auto length = this->getSoundLength();
    while (tracker_l < length) {
        // Get frame counter
        left_frame = tracker_l * _sndinfo.samplerate;
        right_frame = (tracker_l + beat_length) * (double) _sndinfo.samplerate;

        swingFrames(left_frame, right_frame);

        // Update tracker
        tracker_l = offset + beat_length * ++beat;
    }
}

inline double SoundFile::getSoundLength() const {
    return (double) _sndinfo.frames / _sndinfo.samplerate;
}

void SoundFile::swingFrames(const uint32_t leftFrame,
                            const uint32_t rightFrame) {
    return;
}
