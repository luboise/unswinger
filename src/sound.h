#pragma once

#include <sndfile.h>
#include <string>
#include <vector>

class SoundFile {
   public:
    SoundFile(const std::string& filepath);
    ~SoundFile();

    void addSwing(const double bpm, double offset);

    inline bool isValid() const;
    inline double getSoundLength() const;

   private:
    std::vector<double> _samples;
    SF_INFO _sndinfo;

    void swingFrames(const uint32_t leftFrame, const uint32_t rightFrame);
};