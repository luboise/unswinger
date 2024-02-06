#include "sound.h"

#include <complex>
#include <filesystem>
#include <iostream>

#define SWING_RATIO (2.0 / 3.0)
#define PITCH_WINDOW (44100 / 10)

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

SampleList SoundFile::getStretched(double stretchFactor) const {
    return getStretched(stretchFactor, 0, _sndinfo.frames - 1);
}

SampleList SoundFile::getStretched(double stretchFactor, int32_t startFrame,
                                   int32_t endFrame) const {
    const int32_t startSample = startFrame * _sndinfo.channels;

    const int32_t frame_count = endFrame - startFrame;
    const int32_t sample_count = frame_count * _sndinfo.channels;

    int32_t out_frames = floor(frame_count / stretchFactor);

    SampleList stretched(out_frames * _sndinfo.channels);

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
            double interpolated_val =
                weight_l * left_frame[i] + weight_r * right_frame[i];
            stretched[current_sample++] = interpolated_val;
        }
    }

    return stretched;
}

std::vector<SampleType> SoundFile::getFrame(const int32_t frameIndex) const {
    std::vector<SampleType> frame(_sndinfo.channels);

    auto startSample = frameIndex * _sndinfo.channels;
    for (size_t i = 0; i < _sndinfo.channels; i++) {
        frame[i] = _samples[startSample + i];
    }

    return frame;
}
void SoundFile::setSamples(const SampleList samples) {
    this->_samples = samples;
    this->_sndinfo.frames = ceilf((double)samples.size() / _sndinfo.channels);
}
void SoundFile::setChannel(const size_t channel, const size_t offset,
                           const SampleList samples) {
    size_t index;
    for (size_t i = 0; i < samples.size(); i++) {
        index = (offset + i) * _sndinfo.channels + channel;
        _samples[index] = samples[i];
    }
}

SampleList SoundFile::getSamples() const { return this->_samples; }

SampleList SoundFile::getChannel(uint8_t channel) const {
    SampleList channel_data(_sndinfo.frames);
    for (size_t i = 0; i < _sndinfo.frames; i++) {
        channel_data[i] = _samples[i * _sndinfo.channels + channel];
    }

    return channel_data;
}

SampleList SoundFile::getPitched(const SampleList& channelData,
                                 const double semitones) const {
    SampleList outputList(channelData.size());

    for (size_t i = 0; i < channelData.size(); i += PITCH_WINDOW) {
        auto chunk = SampleList(
            channelData.begin() + 1 + i,
            channelData.begin() +
                std::clamp(i + PITCH_WINDOW, (size_t)0, channelData.size()));
        changePitch(chunk, semitones);
        for (size_t j = 0; j < chunk.size(); j++) {
            outputList[i + j] = chunk[j];
        }
    }

    return outputList;
}

void SoundFile::changePitch(SampleList& inplaceData,
                            const double& semitones) const {
    return changePitch(inplaceData, semitones, 0, inplaceData.size() - 1);
}

void SoundFile::changePitch(SampleList& inplaceData, const double& semitones,
                            const size_t startOffset,
                            const size_t endOffset) const {
    if (startOffset > endOffset) {
        throw std::range_error("Start offset is after end offset.");
    }

    if (startOffset >= inplaceData.size() || endOffset >= inplaceData.size()) {
        throw std::domain_error("Offsets are out of range of sample list.");
    }

    SampleList channelBackup(inplaceData.begin() + 1 + startOffset,
                             inplaceData.begin() + 1 + endOffset);

    std::vector<std::complex<SampleType>> complexData(channelBackup.size());

    fftw_plan plan = fftw_plan_dft_r2c_1d(
        complexData.size(), channelBackup.data(),
        reinterpret_cast<fftw_complex*>(complexData.data()), FFTW_ESTIMATE);
    fftw_execute(plan);

    double binSize = (double)_sndinfo.samplerate / channelBackup.size();

    // New array to copy into
    std::vector<std::complex<double>> complexModifiedData(complexData.size());

    // Multiplier = a 12th of 2 (want double to be an octave)
    double semitoneMultiplier = powf(2.0, abs(semitones) / 12.0);
    if (semitones < 0.0) {
        semitoneMultiplier = 1 / semitoneMultiplier;
    }

    uint32_t nyquistFrequency = _sndinfo.samplerate / 2;

    // Shift everything up by wanted amount of semitones
    const size_t size = complexData.size();
    for (size_t i = 1; i < size; i++) {
        double new_index = i * semitoneMultiplier;
        double frequency = new_index * binSize;
        if (frequency >= nyquistFrequency) {
            break;
        }

        if (new_index < size) {
            complexModifiedData[(uint32_t)new_index] = complexData[i];
        }
    }

    fftw_plan c2r_plan = fftw_plan_dft_c2r_1d(
        complexData.size(),
        reinterpret_cast<fftw_complex*>(complexModifiedData.data()),
        channelBackup.data(), FFTW_ESTIMATE);
    fftw_execute(c2r_plan);

    fftw_destroy_plan(plan);
    fftw_destroy_plan(c2r_plan);

    for (auto& value : channelBackup) {
        value /= channelBackup.size();
    }

    for (size_t i = 0; i < endOffset - startOffset; i++) {
        inplaceData[i + startOffset] = channelBackup[i];
    }
}

void SoundFile::addSwingFourier(const double bpm, double offset) {
    // Modify Pitch
    double beat_length;
    beat_length = 60 / bpm;

    if (offset > 0) {
        offset = std::fmod(offset, beat_length) - beat_length;
    }

    for (size_t channel_index = 0; channel_index < _sndinfo.channels;
         channel_index++) {
        SampleList channelData = this->getChannel(channel_index);

        uint32_t beat = 0;

        // Trackers in seconds
        double tracker_l = offset;

        int32_t left_frame = -1;
        int32_t right_frame = -1;

        const auto length = this->getSoundLength();
        while (tracker_l < length) {
            // Get frame counter
            left_frame = tracker_l * _sndinfo.samplerate;
            right_frame =
                (tracker_l + beat_length) * (double)_sndinfo.samplerate;

            makeSwung(channelData, left_frame, right_frame);

            // Update tracker
            tracker_l = offset + beat_length * ++beat;
        }

        this->setChannel(channel_index, 0, channelData);
    }
};

void SoundFile::changeVolume(const double newVolume) {
    for (auto& sample : _samples) {
        sample *= newVolume;
    }
}

void SoundFile::normalise() {
    if (_samples.size() == 0) return;
    double maxSample = abs(_samples[0]);
    for (const auto& sample : _samples) {
        auto newVal = abs(sample);
        if (newVal > maxSample) maxSample = newVal;
    }

    double scalar = 1 / maxSample;
    for (auto& sample : _samples) sample *= scalar;
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

void SoundFile::makeSwung(SampleList& samples, uint32_t leftFrame,
                          uint32_t rightFrame) const {
    const size_t frame_count = rightFrame - leftFrame;
    if (leftFrame > rightFrame) {
        throw std::logic_error("Right frame is above left frame.");
    }

    if (leftFrame >= samples.size() || rightFrame >= samples.size()) {
        return;
        // TODO: ACCOUNT FOR FINAL FRAME AND DONT JUST RETURN
        // throw std::range_error("Left or right frame is out of bounds.");
    }

    auto semitonesUp = 12 * log2(2 * SWING_RATIO);
    auto semitonesDown = 12 * log2(2 * (1 - SWING_RATIO));

    // Perform pitch changes
    size_t end = rightFrame;
    size_t middle = (frame_count / 2) + leftFrame;
    changePitch(samples, semitonesUp, leftFrame, middle);
    changePitch(samples, semitonesDown, middle + 1, rightFrame);

    std::vector<double> splice(frame_count);
    splice[0] = samples[0];
    splice[splice.size() - 1] = samples[splice.size() - 1];

    for (size_t frame_index = 1; frame_index < frame_count - 1; frame_index++) {
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

        int64_t leftSample = transposed_pos * frame_count + leftFrame;

        double val = (weight_l * samples[leftSample]) +
                     (weight_r * samples[leftSample + 1]);

        splice[frame_index] = val;
    }

    for (size_t i = 0; i < splice.size(); i++) {
        samples[leftFrame + i] = splice[i];
    }
}