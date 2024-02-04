#include "sound.h"

#include <iostream>

#include <filesystem>

#define SWING_RATIO (2.0 / 3.0)

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

std::vector<double> SoundFile::getStretched(double stretchFactor) const {
    return getStretched(stretchFactor, 0, _sndinfo.frames - 1);
}

std::vector<double> SoundFile::getStretched(double stretchFactor,
                                            int32_t startFrame,
                                            int32_t endFrame) const {
    const int32_t startSample = startFrame * _sndinfo.channels;

    const int32_t frame_count = endFrame - startFrame;
    const int32_t sample_count = frame_count * _sndinfo.channels;

    int32_t out_frames = floor(frame_count / stretchFactor);

    std::vector<double> stretched(out_frames * _sndinfo.channels);

    uint32_t current_sample = 0;
    for (int frame = 0; frame < out_frames; frame++) {
        // pos is [0.0, 1.0]
        double pos = (double)frame / (frame_count - 1);

        // Get the frame offset that we want the real value from
        double offset = pos * frame_count * stretchFactor;

        double weight_r = fmod(offset, 1.0);
        double weight_l = 1 - weight_r;

        int32_t wanted_frame = (startFrame + (uint32_t)offset);
        auto left_frame = this->getFrame(wanted_frame);
        auto right_frame = this->getFrame(wanted_frame + 1);

        for (size_t i = 0; i < _sndinfo.channels; i++) {
            stretched[current_sample++] =
                weight_l * left_frame[i] + weight_r * right_frame[i];
        }
    }

    return stretched;
}

std::vector<double> SoundFile::getFrame(const int32_t frameIndex) const {
    std::vector<double> frame(_sndinfo.channels);

    auto startSample = frameIndex * _sndinfo.channels;
    for (size_t i = 0; i < _sndinfo.channels; i++) {
        frame[i] = _samples[startSample + i];
    }

    return frame;
}
void SoundFile::setSamples(const std::vector<double> samples) {
    this->_samples = samples;
    this->_sndinfo.frames = ceilf((double)samples.size() / _sndinfo.channels);
};

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

    for (size_t current_channel = 0; current_channel < _sndinfo.channels;
         current_channel++) {
        for (size_t frame_index = 1; frame_index < frame_count - 1;
             frame_index++) {
            double pos = (double)frame_index / (double)(frame_count - 1);

            // Normalised pos
            double transposed_pos;
            if (pos <= SWING_RATIO) {
                transposed_pos = pos / SWING_RATIO * 0.5;
            } else {
                transposed_pos = (pos - SWING_RATIO) / (2 * (1 - SWING_RATIO));
                transposed_pos += 0.5;
            }

            // suppose relevantframe = 3.4, then we want 60% of 3 and 40% of 4
            double weight_r = fmod(transposed_pos, 1.0);
            double weight_l = 1 - weight_r;

            int64_t localFrame = transposed_pos * frame_count;
            int64_t localSample =
                localFrame * _sndinfo.channels + current_channel;

            int64_t lsampleindex = leftSample + localSample;
            int64_t rsampleindex = lsampleindex + _sndinfo.channels;

            if (lsampleindex < 0 || lsampleindex >= global_sample_count ||
                rsampleindex >= global_sample_count) {
                continue;
            }

            double val = (weight_l * _samples[lsampleindex]) +
                         (weight_r * _samples[rsampleindex]);

            auto test = frame_index * _sndinfo.channels + current_channel;
            splice[test] = val;
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
