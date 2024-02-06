#pragma once

#include <sndfile.h>
#include <array>
#include <string>
#include <vector>

#include <fftw3.h>

typedef double SampleType;
typedef std::vector<SampleType> SampleList;

class SoundFile {
   public:
    SoundFile(const std::string& filepath);
    ~SoundFile();

    void addSwing(const double bpm, double offset);
    void addSwingFourier(const double bpm, double offset);

    void exportToFile(const std::string filename) const;

    inline bool isValid() const;
    inline double getSoundLength() const;

    SampleList getStretched(double stretchFactor) const;
    SampleList getStretched(double stretchFactor, int32_t startFrame,
                            int32_t endFrame) const;
    std::vector<SampleType> getFrame(const int32_t frameIndex) const;
    SampleList getSamples() const;
    SampleList getChannel(uint8_t channel) const;

    SampleList getPitched(const std::vector<double>& channelData,
                          const double semitones) const;

    void setSamples(const std::vector<double> frames);
    void setChannel(const size_t channel, const size_t offset,
                    const std::vector<double> samples);

    void changeVolume(const double newVolume);
    void normalise();

   private:
    SampleList _samples;
    SF_INFO _sndinfo;

    void swingFrames(const int32_t leftFrame, const int32_t rightFrame);
    void makeSwung(SampleList& samples, uint32_t leftFrame,
                   uint32_t rightFrame) const;

    void changePitch(SampleList& inplaceData, const double& semitones) const;
    void changePitch(SampleList& inplaceData, const double& semitones,
                     const size_t startOffset, const size_t endOffset) const;
};