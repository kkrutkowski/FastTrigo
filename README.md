FastTrigo 1.0 (c) 2013 Robin Lobel
Modified FMA/AVX512 implementations
=========
Fast yet accurate trigonometric functions

Namespace FTA has 2 sets of functions:

    Packed Scalar: sin/cos functions computing 4 to 16 values at the same time (using FMA/AVX512 instruction sets)
  
Max error:

    FTA::sin_ps/cos_ps/sin_2pi_ps/cos_2pi_ps - < 1e-5
    FTA::sin_pd/cos_pd/sin_2pi_pd/cos_2pi_pd - < 1e-12

FTA Speed up (FMA, AMD Ryzen 4600h, g++ 13.2, -O3):

    FTA::sin_ps/cos_ps max speed up: x25 (from standard sin + cos)
    FTA::sin_2pi_ps/cos_2pi_ps max speed up: x32 (from standard sin + cos)

Distributed under Revised BSD License
