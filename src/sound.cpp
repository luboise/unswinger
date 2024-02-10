#include "sound.h"
#include "utils.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <algorithm>
#include <complex>
#include <filesystem>
#include <iostream>

#define SWING_RATIO (2.0 / 3.0)
#define PITCH_WINDOW (4096)
#define LANCZOS_WINDOW 3

#define WINDOW_SIZE (4096)
#define OVERLAP_RATIO_FAST (0.25)
#define OVERLAP_RATIO_SLOW (0.6)

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
        SampleList samples(sample_count);

        auto read_count = sf_readf_double(snd, samples.data(), info.frames);
        if (read_count != info.frames) {
            std::cerr << "Read failed for file \"" << filepath
                      << "\", frames read: " << read_count << "/" << info.frames
                      << std::endl;
            return;
        }

        std::vector<SampleList> channels(info.channels);

        for (size_t i = 0; i < samples.size(); i++) {
            channels[i % info.channels].push_back(samples[i]);
        }

        this->_sndinfo = info;
        this->_channels = channels;

        sf_close(snd);

    } else {
        std::cerr << "Error opening sound file \"" << filepath << "\""
                  << std::endl;
    }
}

SoundFile::~SoundFile() {}

inline bool SoundFile::isValid() const { return this->_channels.size() > 0; }

// void SoundFile::addSwing(const double bpm, double offset) {
//     double beat_length;
//
//     beat_length = 60 / bpm;
//
//     if (offset > 0) {
//         offset = std::fmod(offset, beat_length) - beat_length;
//     }
//
//     uint32_t beat = 0;
//
//     // Trackers in seconds
//     double tracker_l = offset;
//
//     int32_t left_frame = -1;
//     int32_t right_frame = -1;
//
//     const auto length = this->getSoundLength();
//     while (tracker_l < length) {
//         // Get frame counter
//         left_frame = tracker_l * _sndinfo.samplerate;
//         right_frame = (tracker_l + beat_length) *
//         (double)_sndinfo.samplerate;
//
//         swingFrames(left_frame, right_frame);
//
//         // Update tracker
//         tracker_l = offset + beat_length * ++beat;
//     }
// }

void SoundFile::exportToFile(const std::string filename) {
    if (!this->isValid()) {
        std::cerr << "Attempting to write invalid file. Skipping..."
                  << std::endl;
        return;
    }

    this->resizeChannels();

    SampleList outSamples(this->getSampleCount());

    size_t channel_count = _channels.size();

    for (size_t i = 0; i < this->_channels[0].size(); i++) {
        for (size_t j = 0; j < channel_count; j++) {
            outSamples[i * channel_count + j] = _channels[j][i];
        }
    }

    // Remove const binding from const function
    SF_INFO info = _sndinfo;
    SNDFILE* snd = sf_open(filename.c_str(), SFM_WRITE, &info);

    sf_write_double(snd, outSamples.data(), outSamples.size());
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

SampleList SoundFile::getFrame(const int32_t frameIndex) const {
    SampleList frame(_sndinfo.channels);

    for (size_t i = 0; i < _sndinfo.channels; i++) {
        frame[i] = _channels[i][frameIndex];
    }

    return frame;
}

void SoundFile::setChannel(const size_t channel, const size_t offset,
                           const SampleList samples) {
    _channels[channel] = samples;

    /*if (samples.size() <= _channels[channel].size())
    size_t index;
    for (size_t i = 0; i < samples.size(); i++) {
        index = (offset + i) * _sndinfo.channels + channel;
        _channels[channel] = samples[i];
    }*/
}

// SampleList SoundFile::getSamples() const { return this->_samples; }

SampleList SoundFile::getChannel(uint8_t channel) const {
    if (channel >= _channels.size()) {
        std::cerr << "Attempted to get channel at invalid index " << channel
                  << "." << std::endl;
        return SampleList();
    }

    return _channels[channel];
}

SampleList SoundFile::getPitched(const SampleList& channelData,
                                 const double semitones) const {
    SampleList outputList(channelData.size());

    for (size_t i = 0; i < channelData.size(); i += PITCH_WINDOW) {
        auto chunk = SampleList(
            channelData.begin() + i,
            channelData.begin() +
                std::clamp(i + PITCH_WINDOW, (size_t)0, channelData.size()));
        changePitch(chunk, semitones);
        for (size_t j = 0; j < chunk.size(); j++) {
            outputList[i + j] = chunk[j];
        }
    }

    return outputList;
}

size_t SoundFile::getChannelCount() const { return _sndinfo.channels; }
size_t SoundFile::getSampleCount() const {
    if (_channels.size() == 0) return 0;
    return _channels.size() * _channels[0].size();
    ;
}

void SoundFile::changePitch(SampleList& inplaceData,
                            const double& semitones) const {
    return changePitch(inplaceData, semitones, 0, inplaceData.size() - 1);
}

double lanczosWindow(double x, double a) {
    if (std::abs(x) < 1e-8) return 1.0;
    if (std::abs(x) >= a) return 0.0;
    return std::sin(M_PI * x) * std::sin(M_PI * x / a) / (M_PI * M_PI * x * x);
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

    SampleList channelSlice(inplaceData.begin() + startOffset,
                            inplaceData.begin() + 1 + endOffset);

    auto complexData = getFFT(channelSlice, _sndinfo.samplerate / 4);
    double binSize = (double)_sndinfo.samplerate / channelSlice.size();

    // Multiplier = a 12th of 2 (want double to be an octave)
    double semitoneMultiplier = powf(2.0, semitones / 12.0);

    uint32_t nyquistFrequency = _sndinfo.samplerate / 2;

    // Shift everything up by wanted amount of semitones
    const size_t size = complexData.size();

    std::complex<double> complexSum;

    for (size_t i = 0; i < size / 2; i++) {
        complexSum = 0;
        double old_frequency = i * binSize;
        double old_index = i / semitoneMultiplier;
        double new_index = i * semitoneMultiplier;

        double frequency = i * binSize;
        if (frequency >= nyquistFrequency) {
            break;
        }

        double highRatio = fmod(i, 1.0);
        if (i < size) {
            complexData[(uint32_t)new_index] +=
                (1.0 - highRatio) * complexData[i];
            complexData[(uint32_t)new_index + 1] += (highRatio)*complexData[i];
        }

        // for (size_t current_bin = (int32_t)i - LANCZOS_WINDOW + 1;
        //      current_bin <= i + LANCZOS_WINDOW; ++current_bin) {
        //     if (current_bin < 0 || current_bin > size / 2) continue;

        //     double sincArg = current_bin - frequency / binSize;
        //     double weight = lanczosWindow(sincArg, LANCZOS_WINDOW);
        //     //*complexData[old_index].real();

        //     complexSum += weight * complexData[current_bin];
        // }
        // complexModifiedData[i] = complexSum;
    }

    SampleList iFFTValues = getIFFT(complexData);

    for (size_t i = 0; i < endOffset - startOffset; i++) {
        inplaceData[i + startOffset] = iFFTValues[i];
    }
}

SampleList SoundFile::getHamming(const SampleList& samples) {
    SampleList hammingValues(samples.size());

    // Magic numbers for deciding tradeoff
    const double ALPHA = 0.54;
    const double BETA = 0.46;

    const size_t SIZE = samples.size();

    for (size_t i = 0; i < SIZE; i++) {
        double hammingValue =
            ALPHA - BETA * std::cos(2 * M_PI * (double)i / (SIZE - 1));
        hammingValues[i] = hammingValue * samples[i];
    }

    return hammingValues;
}

SampleList SoundFile::getVocoded(const SampleList& samples,
                                 const double stretchFactor) const {
    // std::cout << "Beginning vocode." << std::endl;

    double hopA =
        floor(WINDOW_SIZE * (1 - (stretchFactor > 1) ? OVERLAP_RATIO_FAST
                                                     : OVERLAP_RATIO_SLOW));
    double hopS = 1 / stretchFactor * hopA;

    // std::cout << "Creating windows of original signal." << std::endl;
    //  Get windows (each window is N samples)
    std::vector<SampleList> windows = getWindows(samples, WINDOW_SIZE, hopA);
    // std::cout << "Created windows of original signal.\n" << std::endl;

    // std::cout << "Acquiring FFT data for all bins..." << std::endl;
    //  Get bin data for each window
    auto windowBins = std::vector<FFTBinList>(windows.size());
    for (size_t i = 0; i < windows.size(); i++) {
        // std::cout << "Acquiring FFT data for bin " << i+1 << "...";
        const auto& window = windows[i];
        windowBins[i] = getFFT(window);
    }
    // std::cout << "Acquired FFT bin data.\n" << std::endl;

    double deltaTimeAnalysed = hopA / _sndinfo.samplerate;
    double deltaTimeSynthesised = hopS / _sndinfo.samplerate;

    std::vector<double> binPhasesAnalysed(WINDOW_SIZE);
    std::vector<double> binPhasesSynthesised(WINDOW_SIZE);

    // std::cout << "Synthesising frame data..." << std::endl;

    size_t frames = windows.size();
    // For each window
    for (size_t frame_index = 0; frame_index < frames; frame_index++) {
        // std::cout << "Synthesising frame " << frame_index + 1 << "..."
        //<< std::endl;
        SampleList& frame = windows[frame_index];
        FFTBinList& currentBins = windowBins[frame_index];

        size_t currentBinSize = currentBins.size();

        // For each bin within each window
        for (size_t bin_index = 0; bin_index < currentBins.size();
             bin_index++) {
            FFT_T currentBin = currentBins[bin_index];

            double binFrequency = ((double)bin_index * _sndinfo.samplerate) /
                                  (currentBinSize - 1);

            SampleType previousPhase = binPhasesAnalysed[bin_index];
            SampleType currentPhase = std::arg(currentBin);

            // Calculate deviation of current bin between now and last frame and
            // wrap it to [-pi, pi]
            double binFrequencyDeviation =
                (currentPhase - previousPhase) / (deltaTimeAnalysed);

            double binTrueFrequency =
                WrapAngle(binFrequencyDeviation - binFrequency) + binFrequency;

            double synthesisedPhase = binPhasesSynthesised[bin_index] +
                                      (deltaTimeSynthesised * binTrueFrequency);

            double binMagnitude = std::abs(currentBin);

            FFT_T synthesisedVal = std::polar(binMagnitude, synthesisedPhase);
            currentBin = synthesisedVal;

            // Cleanup
            binPhasesAnalysed[bin_index] = currentPhase;
            binPhasesSynthesised[bin_index] = synthesisedPhase;
        }

        windows[frame_index] = getHamming(getIFFT(currentBins));
    }

    // std::cout << "Synthesised frame data.\n" << std::endl;

    // std::cout << "Reconstructing vocoded signal..." << std::endl;

    SampleList returnSamples(ceil(samples.size() / stretchFactor));
    // For each window
    for (size_t window_index = 0; window_index < windows.size();
         window_index++) {
        auto& window = windows[window_index];
        // Get window starting position in real array
        uint32_t global_window_start = floor(window_index * hopS);

        // Sum the window values to the return vector
        for (size_t sample_index = 0; sample_index < window.size();
             sample_index++) {
            int32_t current_index = global_window_start + sample_index;
            if (current_index < 0 || current_index >= returnSamples.size())
                continue;

            returnSamples[current_index] += window[sample_index];
        }
    }
    // std::cout << "Reconstructed vocoded signal. Vocoding finished.\n"
    //<< std::endl;

    return returnSamples;
}

std::vector<SampleList> SoundFile::getWindows(const SampleList& samples,
                                              const size_t windowSize,
                                              const double hopSize) const {
    std::vector<SampleList> windows;
    uint32_t l_index;
    uint32_t r_index;

    size_t window_index = 0;

    const size_t size = samples.size();
    while (true) {
        l_index = window_index++ * hopSize;
        r_index = l_index + windowSize;

        if (r_index >= size) {
            SampleList window(samples.begin() + l_index, samples.end());
            windows.push_back(window);
            break;
        } else {
            SampleList window(samples.begin() + l_index,
                              samples.begin() + r_index);
            windows.push_back(window);
        }
    }

    return windows;
}

FFTBinList SoundFile::getFFT(const SampleList& samples) const {
    // Call with no minimum
    return getFFT(samples, 0);
}

FFTBinList SoundFile::getFFT(const SampleList& samples,
                             size_t minimumSize) const {
    size_t FFTSize = std::max(samples.size(), minimumSize);
    FFTBinList complexData(FFTSize);

    SampleList samplesCopy = samples;

    fftw_plan plan = fftw_plan_dft_r2c_1d(
        samples.size(), samplesCopy.data(),
        reinterpret_cast<fftw_complex*>(complexData.data()), FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    return complexData;
}

SampleList SoundFile::getIFFT(const FFTBinList& complexData) const {
    SampleList samples(complexData.size());
    auto complexDataCopy = complexData;

    fftw_plan c2r_plan = fftw_plan_dft_c2r_1d(
        complexData.size(),
        reinterpret_cast<fftw_complex*>(complexDataCopy.data()), samples.data(),
        FFTW_ESTIMATE);
    fftw_execute(c2r_plan);
    fftw_destroy_plan(c2r_plan);

    for (auto& value : samples) {
        value /= samples.size();
    }

    return samples;
}

double SoundFile::getNewBeatPosition(double beat, double ratio) {
    double length_change = 2.0 * ratio;

    double subbeat = fmod(beat, 1.0);

    if (subbeat < ratio) {
        subbeat *= 2 * ratio;
    } else {
        subbeat -= ratio;
        subbeat *= (1 - ratio);
        subbeat += ratio;
    }

    return floor(beat) + subbeat;
}

void SoundFile::resizeChannels() {
    uint32_t maxSize = 0;
    for (const auto& channel : _channels) {
        if (channel.size() > maxSize) maxSize = channel.size();
    }

    for (auto& channel : _channels) {
        if (channel.size() < maxSize) {
            channel.resize(maxSize);
        }
    }

    _sndinfo.frames = (_channels.size() > 0) ? _channels[0].size() : 0;
}

void SoundFile::addSwingFourier(const double bpm, double offset) {
    addSwingFourier(bpm, offset, false);
}

void SoundFile::addSwingFourier(const double bpm, double offset,
                                const bool removeSwing) {
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

            makeSwung(channelData, left_frame, right_frame, removeSwing);

            // Update tracker
            tracker_l = offset + beat_length * ++beat;
        }

        this->setChannel(channel_index, 0, channelData);
    }
};

void SoundFile::addSwingVocoded(const double bpm, double offset) {
    addSwingVocoded(bpm, offset, false);
}

void SoundFile::addSwingVocoded(const double bpm, double offset,
                                const bool removeSwing) {
    // Modify Pitch
    double beat_length;
    beat_length = 60 / bpm;

    if (offset > 0) {
        offset = std::fmod(offset, beat_length) - beat_length;
    }

    double usedRatio = (removeSwing) ? 1 - SWING_RATIO : SWING_RATIO;

    const double upMultiplier = (2 * (1 - usedRatio));
    const double downMultiplier = (2 * usedRatio);

    for (size_t channel_index = 0; channel_index < _sndinfo.channels;
         channel_index++) {
        SampleList channelData = this->getChannel(channel_index);

        SampleList spedData = getVocoded(channelData, 1 / upMultiplier);
        SampleList slowedData = getVocoded(channelData, 1 / downMultiplier);

        SampleList newData(channelData.size());
        for (size_t i = 0; i < newData.size(); i++) {
            // Get the current beat
            double realTime = (double)i / _sndinfo.samplerate;
            realTime += offset;
            double current_beat = realTime / beat_length;

            int32_t beat_count = (int32_t)current_beat;

            double scaledBeat = getNewBeatPosition(current_beat, usedRatio);
            double subbeat = fmod(scaledBeat, 1.0);

            if (subbeat < usedRatio) {
                // Convert beat to index of array
                double getIndex = scaledBeat * 2 * usedRatio * beat_length *
                                  _sndinfo.samplerate;
                if (floor(getIndex + 1) < spedData.size()) {
                    // newData[i] = weight_l * spedData[floor(new_index)] +
                    //              weight_r * spedData[floor(new_index) + 1];
                    newData[i] = spedData[getIndex];
                }

            } else {
                double getIndex = scaledBeat * 2 * (1 - usedRatio) *
                                  beat_length * _sndinfo.samplerate;
                if (floor(getIndex + 1) < slowedData.size()) {
                    // newData[i] = weight_l * slowedData[floor(new_index)] +
                    //              weight_r * slowedData[floor(new_index) + 1];
                    newData[i] = slowedData[getIndex];
                }
            }
        }

        this->setChannel(channel_index, 0, newData);
    }
}

void SoundFile::changeVolume(const double newVolume) {
    for (auto& channel : _channels) {
        for (auto& sample : channel) {
            sample *= newVolume;
        }
    }
}

void SoundFile::normalise() {
    if (_channels.size() == 0 || _channels[0].size() == 0) {
        return;
    }
    double maxSample = abs(_channels[0][0]);
    for (auto& channel : _channels) {
        for (const auto& sample : channel) {
            auto newVal = abs(sample);
            if (newVal > maxSample) maxSample = newVal;
        }
    }

    double scalar = 1 / maxSample;
    for (auto& channel : _channels) {
        for (auto& sample : channel) {
            sample *= scalar;
        }
    }
}

// void SoundFile::swingFrames(const int32_t leftFrame, const int32_t
// rightFrame) {
//     const size_t frame_count = rightFrame - leftFrame;
//     const size_t sample_count = frame_count * _sndinfo.channels;
//     const size_t global_sample_count = _sndinfo.frames *
//     _sndinfo.channels;
//
//     const int32_t leftSample = leftFrame * _sndinfo.channels;
//
//     std::vector<double> splice(sample_count);
//
//     for (size_t i = 0; i < sample_count; i++) {
//         auto index = leftSample + i;
//
//         if (index >= 0 && index < global_sample_count) {
//             splice[i] = _samples[index];
//         } else {
//             splice[i] = 0;
//         }
//     }
//
//     for (size_t current_channel = 0; current_channel < _sndinfo.channels;
//          current_channel++) {
//         for (size_t frame_index = 1; frame_index < frame_count - 1;
//              frame_index++) {
//             double pos = (double)frame_index / (double)(frame_count - 1);
//
//             // Normalised pos
//             double transposed_pos;
//             if (pos <= SWING_RATIO) {
//                 transposed_pos = pos / SWING_RATIO * 0.5;
//             } else {
//                 transposed_pos = (pos - SWING_RATIO) / (2 * (1 -
//                 SWING_RATIO)); transposed_pos += 0.5;
//             }
//
//             // suppose relevantframe = 3.4, then we want 60% of 3 and 40%
//             of 4 double weight_r = fmod(transposed_pos, 1.0); double
//             weight_l = 1 - weight_r;
//
//             int64_t localFrame = transposed_pos * frame_count;
//             int64_t localSample =
//                 localFrame * _sndinfo.channels + current_channel;
//
//             int64_t lsampleindex = leftSample + localSample;
//             int64_t rsampleindex = lsampleindex + _sndinfo.channels;
//
//             if (lsampleindex < 0 || lsampleindex >= global_sample_count
//             ||
//                 rsampleindex >= global_sample_count) {
//                 continue;
//             }
//
//             double val = (weight_l * _samples[lsampleindex]) +
//                          (weight_r * _samples[rsampleindex]);
//
//             auto test = frame_index * _sndinfo.channels +
//             current_channel; splice[test] = val;
//         }
//     }
//
//     int32_t real_index;
//     for (size_t i = 0; i < sample_count; i++) {
//         real_index = leftSample + i;
//         if (real_index < 0)
//             continue;
//         else if (real_index >= global_sample_count)
//             break;
//
//         //_samples[real_index] = (splice[i] != -2.0) ? splice[i] : 0;
//         _samples[real_index] = splice[i];
//     }
//
//     return;
// }

void SoundFile::makeSwung(SampleList& samples, int32_t leftFrame,
                          int32_t rightFrame) const {
    makeSwung(samples, leftFrame, rightFrame, false);
}

void SoundFile::makeSwung(SampleList& samples, int32_t leftFrame,
                          int32_t rightFrame, const bool removeSwing) const {
    const size_t frame_count = rightFrame - leftFrame;
    if (leftFrame > rightFrame) {
        throw std::logic_error("Right frame is above left frame.");
    }

    if (leftFrame >= samples.size() || rightFrame >= samples.size()) {
        return;
        // TODO: ACCOUNT FOR FINAL FRAME AND DONT JUST RETURN
        // throw std::range_error("Left or right frame is out of bounds.");
    }

    if (leftFrame < 0 || rightFrame < 0) {
        // TODO: ACCOUNT FOR FIRST FRAME AND DONT JUST RETURN
        return;
    }

    auto semitonesUp = 12 * log2(2 * SWING_RATIO);
    auto semitonesDown = 12 * log2(2 * (1 - SWING_RATIO));

    // Swap up and down if removing swing
    if (removeSwing) {
        auto temp = semitonesUp;
        semitonesUp = semitonesDown;
        semitonesDown = temp;
    }

    // Perform pitch changes
    size_t end = rightFrame;
    size_t middle = (frame_count / 2) + leftFrame;
    changePitch(samples, semitonesUp, leftFrame, middle);
    changePitch(samples, semitonesDown, middle + 1, rightFrame);

    std::vector<double> splice(frame_count);
    splice[0] = samples[0];
    splice[splice.size() - 1] = samples[splice.size() - 1];

    double usedRatio = (removeSwing) ? 1 - SWING_RATIO : SWING_RATIO;

    for (size_t frame_index = 1; frame_index < frame_count - 1; frame_index++) {
        double pos = (double)frame_index / (double)(frame_count - 1);

        // Normalised pos
        double transposed_pos;
        if (pos <= usedRatio) {
            transposed_pos = pos / usedRatio * 0.5;
        } else {
            transposed_pos = (pos - usedRatio) / (2 * (1 - usedRatio));
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