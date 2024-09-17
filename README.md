FastTrigo 1.0 (c) 2013 Robin Lobel
=========
Fast yet accurate trigonometric functions

Namsepace FTA has 3 sets of functions:

    Scalar: standard trigonometric functions
    Packed Scalar: same functions computing 4 or 8 values at the same time (using SSE4.1/FMA if available)
  

Max error:

    FTA::sin_ps/cos_ps/sin_2pi_ps/cos_2pi_ps - < 1e-5
    FTA::sin/sin_ps max error: 0.0007%

FTA Speed up (FMA, AMD Ryzen 4600h, g++ 13.2, -O3):

    FTA::sin_ps/cos_ps speed max up: x25 (from standard sin + cos)
    FTA::sin_2pi_ps/cos_2pi_ps speed max up: x32 (from standard sin + cos)

Distributed under Revised BSD License
