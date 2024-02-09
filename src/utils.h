#pragma once

#include <complex>

#define _USE_MATH_DEFINES
#include <math.h>

template <typename T>
T WrapAngle(const T val) {
    T returnVal = fmod(val + M_PI, M_2_PI);
    returnVal -= M_PI;

    return returnVal;
}