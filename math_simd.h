#pragma once

//
// NOTE: SIMD
//

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
