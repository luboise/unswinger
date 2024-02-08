#include "utils.h"

std::complex<double> SafeScaleComplex(const std::complex<double>& complex,
                                      const double ratio) {
    double m = sqrtf(ratio * ratio *
                         (complex.real() * complex.real() +
                          complex.imag() * complex.imag()) -
                     (complex.imag() * complex.imag()));

    auto returnVal = complex;
    returnVal._Val[0] *= m;
    return returnVal;
}