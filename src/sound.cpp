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

    int32_t left_frame = -1;
    int32_t right_frame = -1;

    const auto length = this->getSoundLength();
    while (tracker_l < length) {
        // Get frame counter
        left_frame = tracker_l * _sndinfo.samplerate;
        right_frame = (tracker_l + beat_length) * (double)_sndinfo.samplerate;

        swingFrames(left_frame, right_frame);

        // Update tracker
        tracker_l = offset + beat_length * ++beat;
    }
}

void SoundFile::exportToFile(const std::string filename) const {
    if (!this->isValid()) {
        std::cerr << "Attempting to write invalid file. Skipping..."
                  << std::endl;
        return;
    } else if (_samples.size() % _sndinfo.channels != 0 ||
               _samples.size() / _sndinfo.channels != _sndinfo.frames) {
        std::cerr << "Mismatch between number of samples and channel count. "
                     "Unable to write file. Skipping..."
                  << std::endl;
        return;
    }

    // Remove const binding from const function
    SF_INFO info = _sndinfo;
    SNDFILE* snd = sf_open(filename.c_str(), SFM_WRITE, &info);

    sf_write_double(snd, _samples.data(), _samples.size());
    sf_close(snd);
}

inline double SoundFile::getSoundLength() const {
    return (double)_sndinfo.frames / _sndinfo.samplerate;
}

void SoundFile::swingFrames(const int32_t leftFrame, const int32_t rightFrame) {
    const size_t frame_count = rightFrame - leftFrame;
    const size_t sample_count = frame_count * _sndinfo.channels;
    const size_t global_sample_count = _sndinfo.frames * _sndinfo.channels;

    const int32_t leftSample = leftFrame * _sndinfo.channels;

    std::vector<double> splice(sample_count);

    for (size_t i = 0; i < sample_count; i++) {
        auto index = leftSample + i;

        if (index >= 0 && index < global_sample_count) {
            splice[i] = _samples[index];
        } else {
            splice[i] = 0;
        }
    }
    // splice[0] = _samples[0];
    // splice[1] = _samples[1];
    // splice[sample_count - 1] = _samples[sample_count - 1];
    // splice[sample_count - 2] = _samples[sample_count - 2];

    for (size_t current_channel = 0; current_channel < _sndinfo.channels;
         current_channel++) {
        for (size_t frame_index = 1; frame_index < frame_count - 1;
             frame_index++) {
            double pos = (float)frame_index / (frame_count - 1);

            // Normalised pos
            double transposed_pos;
            if (pos <= 0.5) {
                transposed_pos = pos * 3.0 / 4;
            } else {
                transposed_pos = ((3 * pos) - 1) / 2;
            }

            // suppose relevantframe = 3.4, then we want 60% of 3 and 40% of 4
            double weight_r = fmod(transposed_pos, 1.0);
            double weight_l = 1 - weight_r;

            int16_t localFrame = transposed_pos * frame_count;
            int16_t localSample = localFrame * 2 + current_channel;

            auto lsampleindex = leftSample + localSample;
            auto rsampleindex = lsampleindex + _sndinfo.channels;

            if (lsampleindex < 0 || lsampleindex >= global_sample_count ||
                rsampleindex >= global_sample_count) {
                continue;
            }

            double val = (weight_l * _samples[lsampleindex]) +
                         (weight_r * _samples[rsampleindex]);

            splice[localSample] = val;
        }
    }

    int32_t real_index;
    for (size_t i = 0; i < sample_count; i++) {
        real_index = leftSample + i;
        if (real_index < 0)
            continue;
        else if (real_index >= global_sample_count)
            break;

        //_samples[real_index] = (splice[i] != -2.0) ? splice[i] : 0;
        _samples[real_index] = splice[i];
    }

    return;
}
