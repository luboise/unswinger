#pragma once

#include <sndfile.h>
#include <array>
#include <string>
#include <vector>

class SoundFile {
   public:
    SoundFile(const std::string& filepath);
    ~SoundFile();

    void addSwing(const double bpm, double offset);
    void exportToFile(const std::string filename) const;

    inline bool isValid() const;
    inline double getSoundLength() const;

    std::vector<double> getStretched(double stretchFactor) const;
    std::vector<double> getStretched(double stretchFactor, int32_t startFrame,
                                     int32_t endFrame) const;

    std::vector<double> getFrame(const int32_t frameIndex) const;

    void setSamples(const std::vector<double> frames);

   private:
    std::vector<double> _samples;
    SF_INFO _sndinfo;

    void swingFrames(const int32_t leftFrame, const int32_t rightFrame);
};