#pragma once

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

struct v1u_x4
{
    union
    {
        struct
        {
            __m128i x;
        };

        u32 e[4];
    };

    inline u32 operator[](u32 Index);
};

union v1i_x4
{
    struct
    {
        __m128i x;
    };

    i32 e[4];
};

union v1_x4
{
    struct
    {
        __m128 x;
    };

    f32 e[4];
};

union v2_x4
{
    struct
    {
        v1_x4 x, y;
    };

    v1_x4 e[2];
};

union v3_x4
{
    struct
    {
        v1_x4 x, y, z;
    };

    struct
    {
        v2_x4 xy;
        v1_x4 Ignored0_;
    };

    struct
    {
        v1_x4 Ignored1_;
        v2_x4 yz;
    };

    struct
    {
        v1_x4 r, g, b;
    };

    v1_x4 e[3];
};

union v4_x4
{
    struct
    {
        v1_x4 x, y, z, w;
    };

    struct
    {
        v3_x4 xyz;
        v1_x4 Ignored0_;
    };

    struct
    {
        v2_x4 xy;
        v2_x4 zw;
    };

    struct
    {
        v1_x4 Ignored1_;
        v3_x4 yzw;
    };

    struct
    {
        v1_x4 r, g, b, a;
    };

    struct
    {
        v3_x4 rgb;
        v1_x4 Ignored2_;
    };

    v1_x4 e[4];
};

// NOTE: Matrices are stored column order
struct m2_x4
{
    v2_x4 v[2];
};

// NOTE: Matrices are stored column order
struct m3_x4
{
    v3_x4 v[3];
};

// NOTE: Matrices are stored column order
struct m4_x4
{
    v4_x4 v[4];
};

union q4_x4
{
    struct
    {
        v1_x4 x, y, z, w;
    };
    
    struct
    {
        v3_x4 xyz;
        v1_x4 Ignored0_;
    };
};

struct aabb2_x4
{
    v2_x4 Min;
    v2_x4 Max;
};

struct aabb3_x4
{
    v3_x4 Min;
    v3_x4 Max;
};

// NOTE: Everything C++ needs us to predeclare because it sucks
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
