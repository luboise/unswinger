#pragma once

#include <sndfile.h>
#include <array>
#include <string>
#include <vector>

#include <fftw3.h>

class SoundFile {
   public:
    SoundFile(const std::string& filepath);
    ~SoundFile();

    void addSwing(const double bpm, double offset);
    void addSwingFourier(const double bpm, double offset);

    void exportToFile(const std::string filename) const;

    inline bool isValid() const;
    inline double getSoundLength() const;

    std::vector<double> getStretched(double stretchFactor) const;
    std::vector<double> getStretched(double stretchFactor, int32_t startFrame,
                                     int32_t endFrame) const;
    std::vector<double> getFrame(const int32_t frameIndex) const;
    std::vector<double> getSamples() const;
    std::vector<double> getChannel(uint8_t channel) const;

    void setSamples(const std::vector<double> frames);
    void setChannel(const size_t channel, const size_t offset, const std::vector<double> samples);

   private:
    std::vector<double> _samples;
    SF_INFO _sndinfo;

    void swingFrames(const int32_t leftFrame, const int32_t rightFrame);
};