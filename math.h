#pragma once

// TODO: Remove
#include <math.h>
#include "types.h"

#define MATH_X64 1
//#define MATH_ARM 1

//
// NOTE: Constants
//

#define Pi32 3.14159265359f
#define Tau32 6.28318530717958647692f

//
// NOTE: Scalar
//

union v2
{
    struct
    {
        f32 x, y;
    };

    f32 e[2];
};

union v2i
{
    struct
    {
        i32 x, y;
    };

    i32 e[2];
};

union v3
{
    struct
    {
        f32 x, y, z;        
    };
        
    struct
    {
        v2 xy;
        f32 Ignored0_;
    };

    struct
    {
        f32 Ignored1_;
        v2 yz;
    };

    struct
    {
        f32 r, g, b;
    };

    f32 e[3];
};    

union v4
{
    struct
    {
        f32 x, y, z, w;
    };

    struct
    {
        v3 xyz;
        f32 Ignored0_;
    };

    struct
    {
        v2 xy;
        v2 zw;
    };
    
    struct
    {
        f32 Ignored1_;
        v3 yzw;
    };

    struct
    {
        v3 rgb;
        f32 Ignored2_;
    };
        
    struct
    {
        f32 r, g, b, a;
    };

    f32 e[4];
};

union v4u
{
    struct
    {
        u32 x, y, z, w;
    };

    struct
    {
        v3 xyz;
        u32 Ignored0_;
    };

    struct
    {
        v2 xy;
        v2 zw;
    };
    
    struct
    {
        u32 Ignored1_;
        v3 yzw;
    };

    u32 e[4];
};

// NOTE: Planes

union p3
{
    struct
    {
        v3 Normal;
        f32 d;
    };
};

// NOTE: Matrices are stored column order
union m2
{
    struct
    {
        f32 e[4];
    };

    struct
    {
        v2 v[2];
    };
};

// NOTE: Matrices are stored column order
union m3
{
    struct
    {
        f32 e[12];
    };
    
    struct
    {
        v3 v[3];
    };
};

// NOTE: Matrices are stored column order
union m4
{
    struct
    {
        f32 e[16];
    };

    struct
    {
        v4 v[4];
    };
};

union q4
{
    struct
    {
        f32 x, y, z, w;
    };
    
    struct
    {
        v3 xyz;
        f32 Ignored0_;
    };
};

struct aabb2
{
    v2 Min;
    v2 Max;
};

struct aabb2i
{
    v2i Min;
    v2i Max;
};

struct aabb3
{
    v3 Min;
    v3 Max;
};

struct ray_cast
{
    i16 NumGridsToVisit;
    i16 CurrX;
    i16 CurrY;
    i16 IncX;
    i16 IncY;
    f32 t;
    f32 NextVert;
    f32 NextHoriz;
    v2 DeltaRecip;
};

//
// NOTE: SOA
//

union v2_soa
{
    struct
    {
        f32* x;
        f32* y;
    };

    f32* e[2];
};

union v2i_soa
{
    struct
    {
        i32* x;
        i32* y;
    };

    i32* e[2];
};

union v3_soa
{
    struct
    {
        f32* x;
        f32* y;
        f32* z;        
    };
        
    struct
    {
        v2_soa xy;
        f32* Ignored0_;
    };

    struct
    {
        f32* Ignored1_;
        v2_soa yz;
    };

    struct
    {
        f32* r;
        f32* g;
        f32* b;
    };

    f32* e[3];
};    

union v4_soa
{
    struct
    {
        f32* x;
        f32* y;
        f32* z;
        f32* w;
    };

    struct
    {
        v3_soa xyz;
        f32* Ignored0_;
    };

    struct
    {
        v2_soa xy;
        v2_soa zw;
    };
    
    struct
    {
        f32* Ignored1_;
        v3_soa yzw;
    };

    struct
    {
        v3_soa rgb;
        f32* Ignored2_;
    };
        
    struct
    {
        f32* r;
        f32* g;
        f32* b;
        f32* a;
    };

    f32* e[4];
};

union m2_soa
{
    struct
    {
        f32* e[4];
    };

    struct
    {
        v2_soa v[2];
    };
};

union m3_soa
{
    struct
    {
        f32* e[12];
    };
    
    struct
    {
        v3_soa v[3];
    };
};

union m4_soa
{
    struct
    {
        f32* e[16];
    };

    struct
    {
        v4_soa v[4];
    };
};

union q4_soa
{
    struct
    {
        f32* x;
        f32* y;
        f32* z;
        f32* w;
    };
    
    struct
    {
        v3_soa v;
        f32* Ignored0_;
    };
};

struct aabb2_soa
{
    v2_soa Min;
    v2_soa Max;
};

struct aabb2i_soa
{
    v2i_soa Min;
    v2i_soa Max;
};

struct aabb3_soa
{
    v3_soa Min;
    v3_soa Max;
};

//
// NOTE: SIMD
//

#if MATH_X64

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

#elif MATH_ARM

#include <arm_neon.h>

#endif

struct v1u_x4
{
    union
    {
#if MATH_X64
        struct
        {
            __m128i x;
        };
#elif MATH_ARM
        // TODO: Implement this, rn issues like below
        //struct
        //{
        //    uint32x4_t x;
        //};
#endif
        
        u32 e[4];
    };

    inline u32 operator[](u32 Index);
};

union v1i_x4
{
#if MATH_X64
    struct
    {
        __m128i x;
    };
#elif MATH_ARM
    // TODO: Implement correctly, there are issues with instructions and documentation sucks so get back to this. Rn everythign is scalar
    //struct
    //{
    //    int32x4_t x;
    //};
#endif
    
    i32 e[4];
};
 
union v1_x4
{
#if MATH_X64
    struct
    {
        __m128 x;
    };
#elif MATH_ARM
    // TODO: Implement correctly. THere are some instructions not included in header and others that are maybe not there at all. Get back
    // to this. Rn, everything becomes scalar here
    
    //struct
    //{
    //    float32x4_t x;
    //};
#endif
    
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

//
// NOTE: C++ sucks we need to pre declare
//

inline f32 SquareRoot(f32 A);
inline f32 Sin(f32 A);
inline f32 Tan(f32 A);
inline f32 ArcSin(f32 A);
inline f32 Cos(f32 A);
inline f32 ArcCos(f32 A);
inline f32 ArcTan(f32 X, f32 Y);

inline i8 Abs(i8 A);
inline i16 Abs(i16 A);
inline i32 Abs(i32 A);
inline i64 Abs(i64 A);
inline f32 Abs(f32 A);
inline f64 Abs(f64 A);
inline v2 Abs(v2 A);
inline v3 Abs(v3 A);
inline v4 Abs(v4 A);

inline u8 Min(u8 A, u8 B);
inline u16 Min(u16 A, u16 B);
inline u32 Min(u32 A, u32 B);
inline u64 Min(u64 A, u64 B);
inline i8 Min(i8 A, i8 B);
inline i16 Min(i16 A, i16 B);
inline i32 Min(i32 A, i32 B);
inline i64 Min(i64 A, i64 B);
inline f32 Min(f32 A, f32 B);
inline f64 Min(f64 A, f64 B);
inline v2 Min(v2 A, v2 B);
inline v3 Min(v3 A, v3 B);
inline v4 Min(v4 A, v4 B);

inline u8 Max(u8 A, u8 B);
inline u16 Max(u16 A, u16 B);
inline u32 Max(u32 A, u32 B);
inline u64 Max(u64 A, u64 B);
inline i8 Max(i8 A, i8 B);
inline i16 Max(i16 A, i16 B);
inline i32 Max(i32 A, i32 B);
inline i64 Max(i64 A, i64 B);
inline f32 Max(f32 A, f32 B);
inline f64 Max(f64 A, f64 B);
inline v2 Max(v2 A, v2 B);
inline v3 Max(v3 A, v3 B);
inline v4 Max(v4 A, v4 B);

inline u8 FloorU8(f32 Val);
inline u8 FloorU8(f64 Val);
inline u16 FloorU16(f32 Val);
inline u16 FloorU16(f64 Val);
inline u32 FloorU32(f32 Val);
inline u32 FloorU32(f64 Val);
inline u64 FloorU64(f32 Val);
inline u64 FloorU64(f64 Val);
inline i8 FloorI8(f32 Val);
inline i8 FloorI8(f64 Val);
inline i16 FloorI16(f32 Val);
inline i16 FloorI16(f64 Val);
inline i32 FloorI32(f32 Val);
inline i32 FloorI32(f64 Val);
inline i64 FloorI64(f32 Val);
inline i64 FloorI64(f64 Val);
inline f32 FloorF32(f32 Val);
inline f64 FloorF64(f64 Val);
inline v2 FloorV2(v2 A);
inline v3 FloorV3(v3 A);
inline v4 FloorV4(v4 A);

inline v2 Normalize(v2 A);
inline v3 Normalize(v3 A);
inline v4 Normalize(v4 A);
inline q4 Normalize(q4 A);

inline f32 LengthSquared(v2 A);
inline f32 LengthSquared(v3 A);
inline f32 LengthSquared(v4 A);
inline f32 LengthSquared(q4 Q);

inline f32 Length(v2 A);
inline f32 Length(v3 A);
inline f32 Length(v4 A);
inline f32 Length(q4 Q);

inline f32 Dot(v2 A, v2 B);
inline f32 Dot(v3 A, v3 B);
inline f32 Dot(v4 A, v4 B);

inline m2 M2Identity();
inline m3 M3Identity();
inline m4 M4Identity();

inline m2 Transpose(m2 M);
inline m3 Transpose(m3 M);
inline m4 Transpose(m4 M);

inline m4 M4Pos(v3 Pos);

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

#include "math_scalar.cpp"
#include "math_simd.cpp"
