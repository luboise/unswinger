#pragma once

#include <complex>

#define _USE_MATH_DEFINES
#include <math.h>

template <typename T>
T WrapAngle(const T val) {
    T returnVal = fmod(val + M_PI, 2 * M_PI);
    returnVal -= M_PI;

    return returnVal;
}

template <typename T>
T UnwrapAngle(const T val) {
    T v = val;

    while (v > M_PI) v -= 2 * M_PI;
    while (v < -M_PI) v += 2 * M_PI;

    return v;
}