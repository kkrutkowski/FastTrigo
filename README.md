FastTrigo 1.0 (c) 2013 Robin Lobel
=========
Fast yet accurate trigonometric functions

Namsepace FTA has 3 sets of functions:

    Scalar: standard trigonometric functions
    Packed Scalar: same functions computing 4 or 8 values at the same time (using SSE4.1/AVX if available)
  

FTA Accuracy:

    FTA::sqrt/sqrt_ps max error: 0%
    FTA::atan2/atan2_ps max error: 0.0005%
    FTA::cos/cos_ps max error: 0.0007%
    FTA::sin/sin_ps max error: 0.0007%

FTA Speed up (SSE, MSVC2012 x64):

    FTA::sqrt speed up: x1.5 (from standard sqrt)
    FTA::atan2 speed up: x1.7 (from standard atan2)
    FTA::sin/cos speed up: x1.6 (from standard sin/cos)
    FTA::sincos speed up: x1.8 (from standard sin+cos)
    FTA::sqrt_ps speed up: x4.9 (from standard sqrt)
    FTA::atan2_ps speed up: x5.2 (from standard atan2)
    FTA::sin_ps/cos_ps speed up: x4.3 (from standard sin/cos)
    FTA::sincos_ps speed up: x5.2 (from standard sin+cos)

FTA Speed up (AVX, g++ 13.2, -O3):

    FTA::sincos speed up: x0.63 (from standard sin+cos)
    FTA::sin_ps/cos_ps speed up: x11.2 (from standard sin/cos)
    FTA::sincos_ps speed up: x9.9 (from standard sin+cos)

Distributed under Revised BSD License
