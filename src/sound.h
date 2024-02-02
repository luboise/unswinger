#pragma once

#include <sndfile.h>
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

   private:
    std::vector<double> _samples;
    SF_INFO _sndinfo;

    void swingFrames(const int32_t leftFrame, const int32_t rightFrame);
};