#include "fftplan.h"

#include <complex>
#include <iostream>

Plan::~Plan() { this->destroyPlan(); }

size_t Plan::getSize() const { return _details.planSize; }

void Plan::destroyPlan() {
    if (_details.plan != nullptr) {
        fftw_destroy_plan(_details.plan);
    }
    _details.planSize = 0;
}

FFTPlan::FFTPlan(const size_t planSize) { this->resize(planSize); }

FFTPlan::~FFTPlan() {}

void FFTPlan::resize(const size_t newSize) {
    if (newSize == 0) return;
    this->destroyPlan();

    SampleList inExample(newSize);
    FFTBinList outExample(newSize);

    fftw_plan plan = fftw_plan_dft_r2c_1d(
        newSize, reinterpret_cast<double*>(inExample.data()),
        reinterpret_cast<fftw_complex*>(outExample.data()), FFTW_ESTIMATE);

    if (plan != nullptr) {
        _details.plan = plan;
        _details.planSize = newSize;
    }
}

void FFTPlan::execute(SampleList& in, FFTBinList& out) const {
    if (in.size() < _details.planSize) {
        std::cerr
            << "Attempted to execute FFT plan on samplelist that is too short. "
               "Skipping."
            << std::endl;
        return;
    }

    out.resize(_details.planSize);

    fftw_execute_dft_r2c(_details.plan, reinterpret_cast<double*>(in.data()),
                         reinterpret_cast<fftw_complex*>(out.data()));
}

void IFFTPlan::resize(const size_t newSize) {
    if (newSize == 0) {
        return;
    }
    this->destroyPlan();

    FFTBinList inExample(newSize);
    SampleList outExample(newSize);

    fftw_plan plan = fftw_plan_dft_c2r_1d(
        newSize, reinterpret_cast<fftw_complex*>(inExample.data()),
        outExample.data(), FFTW_ESTIMATE);

    if (plan != nullptr) {
        _details.plan = plan;
        _details.planSize = newSize;
    }
}

void IFFTPlan::execute(FFTBinList& in, SampleList& out) const {
    if (in.size() < _details.planSize) {
        std::cerr
            << "Attempted to execute IFFT plan on samplelist that is too short. "
               "Skipping."
            << std::endl;
        return;
    }

    out.resize(_details.planSize);
    fftw_execute_dft_c2r(_details.plan,
                         reinterpret_cast<fftw_complex*>(in.data()),
                         reinterpret_cast<double*>(out.data()));
    for (auto& value : out) {
        value /= out.size();
    }
}

IFFTPlan::IFFTPlan(const size_t planSize) { this->resize(planSize); }

IFFTPlan::~IFFTPlan() {}