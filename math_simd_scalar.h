#pragma once

#if MATH_SIMD_X64

/* NOTE: Headers: https://stackoverflow.com/questions/11228855/header-files-for-x86-simd-intrinsics
<mmintrin.h>  MMX
<xmmintrin.h> SSE
<emmintrin.h> SSE2
<pmmintrin.h> SSE3
<tmmintrin.h> SSSE3
<smmintrin.h> SSE4.1
<nmmintrin.h> SSE4.2
<ammintrin.h> SSE4A
<wmmintrin.h> AES
<immintrin.h> AVX, AVX2, FMA
*/

#include <xmmintrin.h>
// TODO: Remove
#include <smmintrin.h>
#include <nmmintrin.h>

#elif MATH_SIMD_ARM

#include <arm_neon.h>

#endif

struct v1u_x4
{
    union
    {
#if MATH_SIMD_X64
        struct
        {
            __m128i x;
        };
#elif MATH_SIMD_ARM
        struct
        {
            uint32x4_t x;
        };
#endif
        
        u32 e[4];
    };
};

union v1i_x4
{
#if MATH_SIMD_X64
    struct
    {
        __m128i x;
    };
#elif MATH_SIMD_ARM
    struct
    {
        int32x4_t x;
    };
#endif
    
    i32 e[4];
};
 
union v1_x4
{
#if MATH_SIMD_X64
    struct
    {
        __m128 x;
    };
#elif MATH_SIMD_ARM
    struct
    {
        float32x4_t x;
    };
#endif
    
    f32 e[4];
};

#if 0
inline v1_x4 SquareRoot(v1_x4 A);
inline v2_x4 SquareRoot(v2_x4 A);
inline v3_x4 SquareRoot(v3_x4 A);
inline v4_x4 SquareRoot(v4_x4 A);

inline v1_x4 Sin(v1_x4 x);
inline v2_x4 Sin(v2_x4 Angle);
inline v3_x4 Sin(v3_x4 Angle);
inline v4_x4 Sin(v4_x4 Angle);

inline v1_x4 Cos(v1_x4 x);
inline v2_x4 Cos(v2_x4 Angle);
inline v3_x4 Cos(v3_x4 Angle);
inline v4_x4 Cos(v4_x4 Angle);

inline v1_x4 LengthSquared(v2_x4 A);
inline v1_x4 LengthSquared(v3_x4 A);
inline v1_x4 LengthSquared(v4_x4 A);
inline v1_x4 LengthSquared(q4_x4 Q);
#endif
