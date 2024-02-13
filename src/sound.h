#pragma once

#include <sndfile.h>
#include <array>
#include <complex>
#include <string>
#include <vector>

#include <fftw3.h>

#include "types.h"

class SoundFile {
   public:
    struct FFTParams {
        double overlayWhenFast = 0.55;
        double overlayWhenSlow = 0.7;
        uint16_t fftWindowSize = 2048;
    };

    SoundFile(const std::string& filepath);
    ~SoundFile();

    void addSwing(const double bpm, double offset);
    void addSwingFourier(const double bpm, double offset);
    void addSwingFourier(const double bpm, double offset,
                         const bool removeSwing);

    void addSwingVocoded(const double bpm, double offset);
    void addSwingVocoded(const double bpm, double offset,
                         const bool removeSwing);

    void exportToFile(const std::string filename);

    bool isValid() const;
    double getSoundLength() const;

    SampleList getStretched(double stretchFactor) const;
    SampleList getStretched(double stretchFactor, int32_t startFrame,
                            int32_t endFrame) const;
    std::vector<SampleType> getFrame(const int32_t frameIndex) const;
    // SampleList getSamples() const;
    SampleList getChannel(uint8_t channel) const;

    SampleList getPitched(const SampleList& channelData,
                          const double semitones) const;

    size_t getChannelCount() const;
    size_t getSampleCount() const;

    void setSamples(const SampleList frames);
    void setChannel(const size_t channel, const size_t offset,
                    const SampleList samples);

    void changeVolume(const double newVolume);
    void normalise();

    SampleList getVocoded(const SampleList& samples,
                          const double semitones) const;

    void setFftParams(const FFTParams& params);

   private:
    std::vector<SampleList> _channels;
    SF_INFO _sndinfo;

    static FFTParams _fftParams;

    // void swingFrames(const int32_t leftFrame, const int32_t rightFrame);
    // void swingFrames(const int32_t leftFrame, const int32_t rightFrame,
    // bool removeSwing);

    void makeSwung(SampleList& samples, int32_t leftFrame,
                   int32_t rightFrame) const;
    void makeSwung(SampleList& samples, int32_t leftFrame, int32_t rightFrame,
                   const bool removeSwing) const;

    void changePitch(SampleList& inplaceData, const double& semitones) const;
    void changePitch(SampleList& inplaceData, const double& semitones,
                     const size_t startOffset, const size_t endOffset) const;

    static SampleList getHamming(const SampleList& samples);

    std::vector<SampleList> getWindows(const SampleList& samples,
                                       const size_t windowSize,
                                       const double hopSize) const;

    std::vector<FFTBinList> getWindowBins(
        std::vector<SampleList> windows) const;

    static double getNewBeatPosition(double beat, double ratio,
                                     bool removingSwing);

    void resizeChannels();
};