#ifndef MINTRIG_HPP
#define MINTRIG_HPP

#include <immintrin.h>

//Most ninja tricks used here:
//http://fastcpp.blogspot.fr/2011/03/changing-sign-of-float-values-using-sse.html
//http://www.songho.ca/misc/sse/sse.html
//http://markplusplus.wordpress.com/2007/03/14/fast-sse-select-operation/
//http://www.masmforum.com/board/index.php?PHPSESSID=786dd40408172108b65a5a36b09c88c0&topic=9515.0
//http://cbloomrants.blogspot.fr/2010/11/11-20-10-function-approximation-by_20.html
//http://assemblyrequired.crashworks.org/2009/10/16/timing-square-root/
//http://nghiaho.com/?p=997
//http://www.researchgate.net/publication/3321724_Efficient_approximations_for_the_arctangent_function
//http://www.ganssle.com/approx/approx.pdf
//http://forum.allaboutcircuits.com/newsgroups/viewtopic.php?t=68185

const float invtwopi=0.1591549f;
const float twopi=6.283185f;
const float threehalfpi=4.7123889f;
const float pi=3.141593f;
const float halfpi=1.570796f;
const float quarterpi=0.7853982f;

static const __m256 AVX_SIGNMASK_PS =  _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000));

//SCALAR
namespace FTA{
    __m256 sqrt(__m256 squared);
    __m256 cos_poly(__m256 x);
    __m256 cos(__m256 angle);
    __m256 sin(__m256 angle);
    void sincos(__m256, __m256*, __m256*);

    __m256 sqrt_ps(__m256 squared);
    __m256 cos_poly_ps(__m256 x);
    __m256 cos_ps(__m256 angle);
    __m256 sin_ps(__m256 angle);
    void sincos_ps(__m256, __m256*, __m256*);

    __m256 sqrt_2pi_ps(__m256 squared);
    __m256 cos_2pi_poly_ps(__m256 x);
    __m256 cos_2pi_ps(__m256 angle);
    __m256 sin_2pi_ps(__m256 angle);
    void sincos_2pi_ps(__m256, __m256*, __m256*);
}



__m256 FTA::sqrt_ps(__m256 squared){
    return _mm256_sqrt_ps(squared);
}



// FMA single precision implementation of sin(x) and cos(x)
__m256 FTA::cos_poly_ps(__m256 x) {
    const __m256 c1 = _mm256_set1_ps(9.99996748424514608493904342371928799e-01f);
    const __m256 c2 = _mm256_set1_ps(-4.99926885961409551350860301487563489e-01f);
    const __m256 c3 = _mm256_set1_ps(4.15007402641074982896429010192022451e-02f);
    const __m256 c4 = _mm256_set1_ps(-1.27439641270480759760246083831391933e-03);
    __m256 x2 = _mm256_mul_ps(x, x);

    // Using FMA instructions for more efficient computation
    __m256 result = _mm256_fmadd_ps(c4, x2, c3); // c3 + c4 * x2
    result = _mm256_fmadd_ps(result, x2, c2);    // c2 + (c3 + c4 * x2) * x2
    result = _mm256_fmadd_ps(result, x2, c1);    // c1 + (c2 + (c3 + c4 * x2) * x2) * x2

    return result;
}

__m256 FTA::cos_ps(__m256 angle){
    //clamp to the range 0..2pi

    //fmod(angle,twopi)
    angle=_mm256_sub_ps(angle,_mm256_mul_ps(_mm256_floor_ps(_mm256_mul_ps(angle,_mm256_set1_ps(invtwopi))),_mm256_set1_ps(twopi)));
    //angle = _mm256_fnmadd_ps(_mm256_floor_ps(_mm256_mul_ps(angle, _mm256_set1_ps(invtwopi))), _mm256_set1_ps(twopi), angle); //seemingly slower than AVX

    __m256 cosangle=angle;
    cosangle=_mm256_xor_ps(cosangle, _mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(halfpi), _CMP_GE_OQ), _mm256_xor_ps(cosangle,_mm256_sub_ps(_mm256_set1_ps(pi),angle))));
    cosangle=_mm256_xor_ps(cosangle,_mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(pi), _CMP_GE_OQ), AVX_SIGNMASK_PS));
    cosangle=_mm256_xor_ps(cosangle, _mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(threehalfpi), _CMP_GE_OQ), _mm256_xor_ps(cosangle,_mm256_sub_ps(_mm256_set1_ps(twopi),angle))));

    __m256 result=FTA::cos_poly_ps(cosangle);

    result=_mm256_xor_ps(result,_mm256_and_ps(_mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(halfpi), _CMP_GE_OQ),_mm256_cmp_ps(angle,_mm256_set1_ps(threehalfpi), _CMP_LT_OQ)), AVX_SIGNMASK_PS));
    return result;
}

__m256 FTA::sin_ps(__m256 angle){
    return FTA::cos_ps(_mm256_sub_ps(_mm256_set1_ps(halfpi),angle));
}


void FTA::sincos_ps(__m256 angle, __m256 *sin, __m256 *cos){
    __m256 anglesign=_mm256_or_ps(_mm256_set1_ps(1.f),_mm256_and_ps(AVX_SIGNMASK_PS,angle));
    //clamp to the range 0..2pi

    //fmod(angle,twopi)
    angle=_mm256_sub_ps(angle,_mm256_mul_ps(_mm256_floor_ps(_mm256_mul_ps(angle,_mm256_set1_ps(invtwopi))),_mm256_set1_ps(twopi)));
    // angle = _mm256_fnmadd_ps(_mm256_floor_ps(_mm256_mul_ps(angle, _mm256_set1_ps(invtwopi))), _mm256_set1_ps(twopi), angle); //seemingly slower than AVX


    __m256 cosangle=angle;
    cosangle=_mm256_xor_ps(cosangle, _mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(halfpi), _CMP_GE_OQ), _mm256_xor_ps(cosangle,_mm256_sub_ps(_mm256_set1_ps(pi),angle))));
    cosangle=_mm256_xor_ps(cosangle,_mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(pi), _CMP_GE_OQ),AVX_SIGNMASK_PS));
    cosangle=_mm256_xor_ps(cosangle, _mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(threehalfpi), _CMP_GE_OQ), _mm256_xor_ps(cosangle,_mm256_sub_ps(_mm256_set1_ps(twopi),angle))));

    __m256 result=FTA::cos_poly_ps(cosangle);

    result=_mm256_xor_ps(result,_mm256_and_ps(_mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(halfpi), _CMP_GE_OQ),_mm256_cmp_ps(angle,_mm256_set1_ps(threehalfpi), _CMP_LT_OQ)), AVX_SIGNMASK_PS));
    *cos=result;

    __m256 sinmultiplier=_mm256_mul_ps(anglesign,_mm256_or_ps(_mm256_set1_ps(1.f),_mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(pi), _CMP_GE_OQ),AVX_SIGNMASK_PS)));
    *sin=_mm256_mul_ps(sinmultiplier,FTA::sqrt_ps(_mm256_fnmadd_ps(result, result, _mm256_set1_ps(1.f))));

    return;
}




// FMA single precision implementation of sin(2πx) and cos(2πx)
__m256 FTA::cos_2pi_poly_ps(__m256 x) {
    const __m256 c1 = _mm256_set1_ps(9.99996748424514608493916615939939119e-01f);
    const __m256 c2 = _mm256_set1_ps(-1.97363223756305008622054184241628274e+01f);
    const __m256 c3 = _mm256_set1_ps(6.46807901818389923284666266409658356e+01f);
    const __m256 c4 = _mm256_set1_ps(-7.84122201283542750695287203100567467e+01f);
    __m256 x2 = _mm256_mul_ps(x, x);

    // Using FMA instructions for more efficient computation
    __m256 result = _mm256_fmadd_ps(c4, x2, c3); // c3 + c4 * x2
    result = _mm256_fmadd_ps(result, x2, c2);    // c2 + (c3 + c4 * x2) * x2
    result = _mm256_fmadd_ps(result, x2, c1);    // c1 + (c2 + (c3 + c4 * x2) * x2) * x2

    return result;
}

__m256 FTA::cos_2pi_ps(__m256 angle){
    //clamp to the range 0..1

    //fmod(angle, 1)
    angle = _mm256_sub_ps(angle, _mm256_mul_ps(_mm256_floor_ps(angle), _mm256_set1_ps(1.0f)));

    __m256 cosangle=angle;
    cosangle=_mm256_xor_ps(cosangle, _mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(0.25f), _CMP_GE_OQ), _mm256_xor_ps(cosangle,_mm256_sub_ps(_mm256_set1_ps(0.5),angle))));
    cosangle=_mm256_xor_ps(cosangle,_mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(0.5f), _CMP_GE_OQ), AVX_SIGNMASK_PS));
    cosangle=_mm256_xor_ps(cosangle, _mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(0.75f), _CMP_GE_OQ), _mm256_xor_ps(cosangle,_mm256_sub_ps(_mm256_set1_ps(1),angle))));

    __m256 result=FTA::cos_2pi_poly_ps(cosangle);

    result=_mm256_xor_ps(result,_mm256_and_ps(_mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(0.25f), _CMP_GE_OQ),_mm256_cmp_ps(angle,_mm256_set1_ps(0.75f), _CMP_LT_OQ)), AVX_SIGNMASK_PS));
    return result;
}

__m256 FTA::sin_2pi_ps(__m256 angle){
    return FTA::cos_2pi_ps(_mm256_sub_ps(_mm256_set1_ps(0.25f),angle));
}


void FTA::sincos_2pi_ps(__m256 angle, __m256 *sin, __m256 *cos){
    __m256 anglesign=_mm256_or_ps(_mm256_set1_ps(1.f),_mm256_and_ps(AVX_SIGNMASK_PS, angle));
    //clamp to the range 0..1

    //fmod(angle, 1)
    angle = _mm256_sub_ps(angle, _mm256_mul_ps(_mm256_floor_ps(angle), _mm256_set1_ps(1.0f)));
    // angle = _mm256_fnmadd_ps(_mm256_floor_ps(_mm256_mul_ps(angle, _mm256_set1_ps(invtwopi))), _mm256_set1_ps(twopi), angle); //seemingly slower than AVX


    __m256 cosangle=angle;
    cosangle=_mm256_xor_ps(cosangle, _mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(0.25f), _CMP_GE_OQ), _mm256_xor_ps(cosangle,_mm256_sub_ps(_mm256_set1_ps(0.5f),angle))));
    cosangle=_mm256_xor_ps(cosangle,_mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(0.5f), _CMP_GE_OQ),AVX_SIGNMASK_PS));
    cosangle=_mm256_xor_ps(cosangle, _mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(0.75f), _CMP_GE_OQ), _mm256_xor_ps(cosangle,_mm256_sub_ps(_mm256_set1_ps(1.0f),angle))));

    __m256 result=FTA::cos_poly_ps(cosangle);

    result=_mm256_xor_ps(result,_mm256_and_ps(_mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(0.25f), _CMP_GE_OQ),_mm256_cmp_ps(angle,_mm256_set1_ps(0.75f), _CMP_LT_OQ)), AVX_SIGNMASK_PS));
    *cos=result;

    __m256 sinmultiplier=_mm256_mul_ps(anglesign,_mm256_or_ps(_mm256_set1_ps(1.f),_mm256_and_ps(_mm256_cmp_ps(angle,_mm256_set1_ps(0.5f), _CMP_GE_OQ),AVX_SIGNMASK_PS)));
    *sin=_mm256_mul_ps(sinmultiplier,FTA::sqrt_ps(_mm256_fnmadd_ps(result, result, _mm256_set1_ps(1.f))));

    return;
}


#endif // MINTRIG_HPP
