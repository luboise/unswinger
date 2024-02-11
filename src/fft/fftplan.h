#pragma once

#include "types.h"

#include <fftw3.h>

struct PlanDetails {
    fftw_plan plan;
    size_t planSize;
};

// Abstract class plan
class Plan {
   public:
    virtual ~Plan();

    virtual void resize(const size_t newSize) = 0;

    size_t getSize() const;

   protected:
    void destroyPlan();

    // Stores new plan size to be used by initialise
    PlanDetails _details{};
};

class FFTPlan : public Plan {
   public:
    FFTPlan(const size_t planSize);

    template <typename T>
    FFTPlan(const std::vector<T>& samples) : FFTPlan(samples.size()){};

    ~FFTPlan();

    void resize(const size_t newSize) override;
    void execute(SampleList& in, FFTBinList& out) const;
};

class IFFTPlan : public Plan {
   public:
    IFFTPlan(const size_t planSize);

    template <typename T>
    IFFTPlan(const std::vector<T>& samples) : IFFTPlan(samples.size()){};

    IFFTPlan(const FFTPlan& fft) : IFFTPlan(fft.getSize()){};
    ~IFFTPlan();

    void resize(const size_t newSize) override;
    void execute(FFTBinList& in, SampleList& out) const;
};