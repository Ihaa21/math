

// =======================================================================================================================================
// NOTE: v1u_x4 Simd
// =======================================================================================================================================

//
// NOTE: Init
//

inline v1u_x4 V1UX4(u32 X)
{
    v1u_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_set1_epi32(X);
#elif MATH_SIMD_ARM
    Result.x = vld1q_dup_s32(&X);
#elif MATH_SIMD_TEST
    Result.e[0] = X;
    Result.e[1] = X;
    Result.e[2] = X;
    Result.e[3] = X;
#endif

    return Result;
}

inline v1u_x4 V1UX4(u32 X, u32 Y, u32 Z, u32 W)
{
    v1u_x4 Result = {};
    
#if MATH_SIMD_X64
    Result.x = _mm_set_epi32(W, Z, Y, X);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = X;
    Result.e[1] = Y;
    Result.e[2] = Z;
    Result.e[3] = W;    
#endif
    
    return Result;
}

//
// NOTE: Memory Ops
//

// TODO: Aligned Masked Load?

inline v1u_x4 V1UX4LoadAligned(u32* X)
{
    v1u_x4 Result = {};
#if MATH_SIMD_X64
    Result.x = _mm_load_si128((__m128i*)X);
#elif MATH_SIMD_ARM
    Result.x = vld1q_u32((u32*)X);
#elif MATH_SIMD_TEST
    u32* CurrPtr = (u32*)X;
    Result.e[0] = CurrPtr[0];
    Result.e[1] = CurrPtr[1];
    Result.e[2] = CurrPtr[2];
    Result.e[3] = CurrPtr[3];
    
#endif
    return Result;
}

inline v1u_x4 V1UX4LoadUnAligned(u32* X)
{
    v1u_x4 Result = {};
#if MATH_SIMD_X64
    Result.x = _mm_loadu_si128((__m128i*)X);
#elif MATH_SIMD_ARM
    Result.x = vld1q_u32((u32*)X);
#elif MATH_SIMD_TEST
    u32* CurrPtr = (u32*)X;
    Result.e[0] = CurrPtr[0];
    Result.e[1] = CurrPtr[1];
    Result.e[2] = CurrPtr[2];
    Result.e[3] = CurrPtr[3];
    
#endif
    return Result;
}

inline v1u_x4 V1UX4Gather(u32* X, u32* Y, u32* Z, u32* W)
{
    v1u_x4 Result = V1UX4(*X, *Y, *Z, *W);
    return Result;
}

inline v1u_x4 V1UX4Gather(u32* Base, v1u_x4 PtrOffset)
{
    v1u_x4 Result = V1UX4Gather(Base + PtrOffset.e[0], Base + PtrOffset.e[1], Base + PtrOffset.e[2], Base + PtrOffset.e[3]);
    return Result;
}

inline v1u_x4 V1UX4Gather(u32* Base, v1i_x4 PtrOffset)
{
    v1u_x4 Result = V1UX4Gather(Base + PtrOffset.e[0], Base + PtrOffset.e[1], Base + PtrOffset.e[2], Base + PtrOffset.e[3]);
    return Result;
}

inline v1u_x4 V1UX4Gather(u32* X, u32* Y, u32* Z, u32* W, v1u_x4 Mask)
{
    // TODO: We can convert __m128i to a mask, its probably faster...
    v1u_x4 Result = V1UX4(Mask.e[0] ? *X : 0, Mask.e[1] ? *Y : 0, Mask.e[2] ? *Z : 0, Mask.e[3] ? *W : 0);
    return Result;
}

inline v1u_x4 V1UX4Gather(u32* Ptr, v1u_x4 PtrOffset, v1u_x4 Mask)
{
    v1u_x4 Result = V1UX4Gather(Ptr + PtrOffset.e[0], Ptr + PtrOffset.e[1], Ptr + PtrOffset.e[2], Ptr + PtrOffset.e[3], Mask);
    return Result;
}

inline v1u_x4 V1UX4Gather(u32* Ptr, v1i_x4 PtrOffset, v1u_x4 Mask)
{
    v1u_x4 Result = V1UX4Gather(Ptr + PtrOffset.e[0], Ptr + PtrOffset.e[1], Ptr + PtrOffset.e[2], Ptr + PtrOffset.e[3], Mask);
    return Result;
}

inline void StoreAligned(v1u_x4 V, void* Dest)
{
#if MATH_SIMD_X64
    _mm_store_si128((__m128i*)Dest, V.x);
#elif MATH_SIMD_ARM
    vst1q_s32((u32*)Dest, V.x);
#elif MATH_SIMD_TEST
    Dest[0] = V.e[0];
    Dest[1] = V.e[1];
    Dest[2] = V.e[2];
    Dest[3] = V.e[3];    
#endif
}

inline void StoreUnAligned(v1u_x4 V, void* Dest)
{
#if MATH_SIMD_X64
    _mm_storeu_si128((__m128i*)Dest, V.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Dest[0] = V.e[0];
    Dest[1] = V.e[1];
    Dest[2] = V.e[2];
    Dest[3] = V.e[3];    
#endif
}

inline void Scatter(v1u_x4 V, u32* X, u32* Y, u32* Z, u32* W)
{
    *X = V.e[0];
    *Y = V.e[1];
    *Z = V.e[2];
    *W = V.e[3];
}

// =======================================================================================================================================
// NOTE: v1i_x4 Simd
// =======================================================================================================================================

//
// NOTE: Init
//

inline v1i_x4 V1IX4(i32 X)
{
    v1i_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_set1_epi32(X);
#elif MATH_SIMD_ARM
    Result.x = vld1q_dup_s32(&X);
#elif MATH_SIMD_TEST
    Result.e[0] = X;
    Result.e[1] = X;
    Result.e[2] = X;
    Result.e[3] = X;
#endif

    return Result;
}

inline v1i_x4 V1IX4(i32 X, i32 Y, i32 Z, i32 W)
{
    v1i_x4 Result = {};
    
#if MATH_SIMD_X64
    Result.x = _mm_set_epi32(W, Z, Y, X);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = X;
    Result.e[1] = Y;
    Result.e[2] = Z;
    Result.e[3] = W;    
#endif
    
    return Result;
}

inline v1i_x4 V1IX4(v1_x4 V)
{
    v1i_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_cvtps_epi32(V.x);
#elif MATH_SIMD_ARM
    Result.x = vcvtq_s32_f32(V.x);
#elif MATH_SIMD_TEST
    Result.e[0] = i32(V.e[0]);
    Result.e[1] = i32(V.e[1]);
    Result.e[2] = i32(V.e[2]);
    Result.e[3] = i32(V.e[3]);    
#endif
    
    return Result;
}

//
// NOTE: Memory Ops
//

// TODO: Aligned Masked Load?

inline v1i_x4 V1IX4LoadAligned(i32* X)
{
    v1i_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_load_si128((__m128i*)X);
#elif MATH_SIMD_ARM
    Result.x = vld1q_s32((i32*)X);
#elif MATH_SIMD_TEST
    Result.e[0] = X[0];
    Result.e[1] = X[1];
    Result.e[2] = X[2];
    Result.e[3] = X[3];
#endif

    return Result;
}

inline v1i_x4 V1IX4LoadUnAligned(i32* X)
{
    v1i_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_loadu_si128((__m128i*)X);
#elif MATH_SIMD_ARM
    Result.x = vld1q_s32((i32*)X);
#elif MATH_SIMD_TEST
    Result.e[0] = X[0];
    Result.e[1] = X[1];
    Result.e[2] = X[2];
    Result.e[3] = X[3];
#endif

    return Result;
}

inline v1i_x4 V1IX4Gather(i32* X, i32* Y, i32* Z, i32* W)
{
    v1i_x4 Result = V1IX4(*X, *Y, *Z, *W);
    return Result;
}

inline v1i_x4 V1IX4Gather(i32* Base, v1u_x4 PtrOffset)
{
    v1i_x4 Result = V1IX4Gather(Base + PtrOffset.e[0], Base + PtrOffset.e[1], Base + PtrOffset.e[2], Base + PtrOffset.e[3]);
    return Result;
}

inline v1i_x4 V1IX4Gather(i32* Base, v1i_x4 PtrOffset)
{
    v1i_x4 Result = V1IX4Gather(Base + PtrOffset.e[0], Base + PtrOffset.e[1], Base + PtrOffset.e[2], Base + PtrOffset.e[3]);
    return Result;
}

inline v1i_x4 V1IX4Gather(i32* X, i32* Y, i32* Z, i32* W, v1u_x4 Mask)
{
    // TODO: We can convert __m128i to a mask, its probably faster...
    v1i_x4 Result = V1IX4(Mask.e[0] ? *X : 0, Mask.e[1] ? *Y : 0, Mask.e[2] ? *Z : 0, Mask.e[3] ? *W : 0);
    return Result;
}

inline v1i_x4 V1IX4Gather(i32* Ptr, v1u_x4 PtrOffset, v1u_x4 Mask)
{
    v1i_x4 Result = V1IX4Gather(Ptr + PtrOffset.e[0], Ptr + PtrOffset.e[1], Ptr + PtrOffset.e[2], Ptr + PtrOffset.e[3], Mask);
    return Result;
}

inline v1i_x4 V1IX4Gather(i32* Ptr, v1i_x4 PtrOffset, v1u_x4 Mask)
{
    v1i_x4 Result = V1IX4Gather(Ptr + PtrOffset.e[0], Ptr + PtrOffset.e[1], Ptr + PtrOffset.e[2], Ptr + PtrOffset.e[3], Mask);
    return Result;
}

inline void StoreAligned(v1i_x4 V, void* Dest)
{
#if MATH_SIMD_X64
    _mm_store_si128((__m128i*)Dest, V.x);
#elif MATH_SIMD_ARM
    vst1q_s32((i32*)Dest, V.x);
#elif MATH_SIMD_TEST
    Dest[0] = V.e[0];
    Dest[1] = V.e[1];
    Dest[2] = V.e[2];
    Dest[3] = V.e[3];    
#endif
}

inline void StoreUnAligned(v1i_x4 V, void* Dest)
{
#if MATH_SIMD_X64
    _mm_storeu_si128((__m128i*)Dest, V.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Dest[0] = V.e[0];
    Dest[1] = V.e[1];
    Dest[2] = V.e[2];
    Dest[3] = V.e[3];    
#endif
}

inline void Scatter(v1i_x4 V, i32* X, i32* Y, i32* Z, i32* W)
{
    *X = V.e[0];
    *Y = V.e[1];
    *Z = V.e[2];
    *W = V.e[3];
}

// =======================================================================================================================================
// NOTE: v1_x4 Simd
// =======================================================================================================================================

//
// NOTE: Init
//

inline v1_x4 V1X4(f32 X)
{
    v1_x4 Result = {};
    
#if MATH_SIMD_X64
    Result.x = _mm_set1_ps(X);
#elif MATH_SIMD_ARM
    Result.x = vld1q_dup_f32(&X);
#elif MATH_SIMD_TEST
    Result.e[0] = X;
    Result.e[1] = X;
    Result.e[2] = X;
    Result.e[3] = X;
#endif
    
    return Result;
}

inline v1_x4 V1X4(f32 X, f32 Y, f32 Z, f32 W)
{
    v1_x4 Result = {};
    
#if MATH_SIMD_X64
    Result.x = _mm_set_ps(W, Z, Y, X);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = X;
    Result.e[1] = Y;
    Result.e[2] = Z;
    Result.e[3] = W;
#endif
    
    return Result;
}

inline v1_x4 V1X4(v1u_x4 V)
{
    v1_x4 Result = {};

    // TODO: Is this correct?
#if MATH_SIMD_X64
    Result.x = _mm_cvtepi32_ps(V.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = f32(V.e[0]);
    Result.e[1] = f32(V.e[1]);
    Result.e[2] = f32(V.e[2]);
    Result.e[3] = f32(V.e[3]);
#endif

    return Result;
}

inline v1_x4 V1X4(v1i_x4 V)
{
    v1_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_cvtepi32_ps(V.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = f32(V.e[0]);
    Result.e[1] = f32(V.e[1]);
    Result.e[2] = f32(V.e[2]);
    Result.e[3] = f32(V.e[3]);
#endif

    return Result;
}

//
// NOTE: Memory Ops
//

inline v1_x4 V1X4LoadAligned(f32* X)
{
    v1_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_load_ps(X);
#elif MATH_SIMD_ARM
    Result.x = vld1q_f32(X);
#elif MATH_SIMD_TEST
    Result.e[0] = X[0];
    Result.e[1] = X[1];
    Result.e[2] = X[2];
    Result.e[3] = X[3];
#endif

    return Result;
}

inline v1_x4 V1X4LoadUnAligned(f32* X)
{
    v1_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_loadu_ps(X);
#elif MATH_SIMD_ARM
    Result.x = vld1q_f32(X);
#elif MATH_SIMD_TEST
    Result.e[0] = X[0];
    Result.e[1] = X[1];
    Result.e[2] = X[2];
    Result.e[3] = X[3];
#endif

    return Result;
}

inline v1_x4 V1X4Gather(f32* X, f32* Y, f32* Z, f32* W)
{
    v1_x4 Result = V1X4(*X, *Y, *Z, *W);
    return Result;
}

inline v1_x4 V1X4Gather(f32* Base, v1u_x4 PtrOffset)
{
    v1_x4 Result = V1X4Gather(Base + PtrOffset.e[0], Base + PtrOffset.e[1], Base + PtrOffset.e[2], Base + PtrOffset.e[3]);
    return Result;
}

inline v1_x4 V1X4Gather(f32* Base, v1i_x4 PtrOffset)
{
    v1_x4 Result = V1X4Gather(Base + PtrOffset.e[0], Base + PtrOffset.e[1], Base + PtrOffset.e[2], Base + PtrOffset.e[3]);
    return Result;
}

inline v1_x4 V1X4Gather(f32* X, f32* Y, f32* Z, f32* W, v1u_x4 Mask)
{
    // TODO: We can convert __m128i to a mask, its probably faster...
    v1_x4 Result = V1X4(Mask.e[0] ? *X : 0, Mask.e[1] ? *Y : 0, Mask.e[2] ? *Z : 0, Mask.e[3] ? *W : 0);
    return Result;
}

inline v1_x4 V1X4Gather(f32* Ptr, v1u_x4 PtrOffset, v1u_x4 Mask)
{
    v1_x4 Result = V1X4Gather(Ptr + PtrOffset.e[0], Ptr + PtrOffset.e[1], Ptr + PtrOffset.e[2], Ptr + PtrOffset.e[3], Mask);
    return Result;
}

inline v1_x4 V1X4Gather(f32* Ptr, v1i_x4 PtrOffset, v1u_x4 Mask)
{
    v1_x4 Result = V1X4Gather(Ptr + PtrOffset.e[0], Ptr + PtrOffset.e[1], Ptr + PtrOffset.e[2], Ptr + PtrOffset.e[3], Mask);
    return Result;
}

inline void StoreAligned(v1_x4 V, void* Dest)
{
#if MATH_SIMD_X64
    _mm_store_ps((f32*)Dest, V.x);
#elif MATH_SIMD_ARM
    vst1q_f32((f32*)Dest, V.x);
#elif MATH_SIMD_TEST
    Dest[0] = V.e[0];
    Dest[1] = V.e[1];
    Dest[2] = V.e[2];
    Dest[3] = V.e[3];
#endif
}

inline void StoreUnAligned(v1_x4 V, void* Dest)
{
#if MATH_SIMD_X64
    _mm_storeu_ps((f32*)Dest, V.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Dest[0] = V.e[0];
    Dest[1] = V.e[1];
    Dest[2] = V.e[2];
    Dest[3] = V.e[3];
#endif
}

inline void Scatter(v1_x4 V, f32* X, f32* Y, f32* Z, f32* W)
{
    *X = V.e[0];
    *Y = V.e[1];
    *Z = V.e[2];
    *W = V.e[3];
}

// =======================================================================================================================================
// NOTE: Common Math Operators
// =======================================================================================================================================

//
// NOTE: Convert
//

//
// NOTE: Cast
//

inline v1u_x4 V1UX4Cast(v1_x4 V)
{
    v1u_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_castps_si128(V.x);
#elif MATH_SIMD_ARM
    Result.x = vreinterpretq_s32_f32(V.x);
#elif MATH_SIMD_TEST
    Result.e[0] = ReinterpretI32(V.e[0]);
    Result.e[1] = ReinterpretI32(V.e[1]);
    Result.e[2] = ReinterpretI32(V.e[2]);
    Result.e[3] = ReinterpretI32(V.e[3]);    
#endif
    
    return Result;
}

inline v1i_x4 V1IX4Cast(v1_x4 V)
{
    v1i_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_castps_si128(V.x);
#elif MATH_SIMD_ARM
    Result.x = vreinterpretq_s32_f32(V.x);
#elif MATH_SIMD_TEST
    Result.e[0] = ReinterpretI32(V.e[0]);
    Result.e[1] = ReinterpretI32(V.e[1]);
    Result.e[2] = ReinterpretI32(V.e[2]);
    Result.e[3] = ReinterpretI32(V.e[3]);    
#endif
    
    return Result;
}

//
// NOTE: Add
// 

inline v1u_x4 operator+(v1u_x4 A, v1u_x4 B)
{
    v1u_x4 Result = {};
    
#if MATH_SIMD_X64
    Result.x = _mm_add_epi32(A.x, B.x);
#elif MATH_SIMD_ARM
    Result.x = vaddq_s32(A.x, B.x);
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] + B.e[0];
    Result.e[1] = A.e[1] + B.e[1];
    Result.e[2] = A.e[2] + B.e[2];
    Result.e[3] = A.e[3] + B.e[3];
#endif

    return Result;
}

inline v1i_x4 operator+(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};
    
#if MATH_SIMD_X64
    Result.x = _mm_add_epi32(A.x, B.x);
#elif MATH_SIMD_ARM
    Result.x = vaddq_s32(A.x, B.x);
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] + B.e[0];
    Result.e[1] = A.e[1] + B.e[1];
    Result.e[2] = A.e[2] + B.e[2];
    Result.e[3] = A.e[3] + B.e[3];
#endif

    return Result;
}

inline v1_x4 operator+(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};
    
#if MATH_SIMD_X64
    Result.x = _mm_add_ps(A.x, B.x);
#elif MATH_SIMD_ARM
    Result.x = vaddq_f32(A.x, B.x);
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] + B.e[0];
    Result.e[1] = A.e[1] + B.e[1];
    Result.e[2] = A.e[2] + B.e[2];
    Result.e[3] = A.e[3] + B.e[3];
#endif
    
    return Result;
}

inline v1u_x4& operator+=(v1u_x4& A, v1u_x4 B)
{
    A = A + B;
    return A;
}

inline v1i_x4& operator+=(v1i_x4& A, v1i_x4 B)
{
    A = A + B;
    return A;
}

inline v1_x4& operator+=(v1_x4& A, v1_x4 B)
{
    A = A + B;
    return A;
}

//
// NOTE: Sub
// 

inline v1u_x4 operator-(v1u_x4 A, v1u_x4 B)
{
    v1u_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_sub_epi32(A.x, B.x);
#elif MATH_SIMD_ARM
    Result.x = vsubq_s32(A.x, B.x);
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] - B.e[0];
    Result.e[1] = A.e[1] - B.e[1];
    Result.e[2] = A.e[2] - B.e[2];
    Result.e[3] = A.e[3] - B.e[3];
#endif

    return Result;
}

inline v1i_x4 operator-(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_sub_epi32(A.x, B.x);
#elif MATH_SIMD_ARM
    Result.x = vsubq_s32(A.x, B.x);
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] - B.e[0];
    Result.e[1] = A.e[1] - B.e[1];
    Result.e[2] = A.e[2] - B.e[2];
    Result.e[3] = A.e[3] - B.e[3];
#endif

    return Result;
}

inline v1_x4 operator-(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};
    
#if MATH_SIMD_X64
    Result.x = _mm_sub_ps(A.x, B.x);
#elif MATH_SIMD_ARM
    Result.x = vsubq_f32(A.x, B.x);
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] - B.e[0];
    Result.e[1] = A.e[1] - B.e[1];
    Result.e[2] = A.e[2] - B.e[2];
    Result.e[3] = A.e[3] - B.e[3];
#endif
    
    return Result;
}

inline v1u_x4& operator-=(v1u_x4& A, v1u_x4 B)
{
    A = A - B;
    return A;
}

inline v1i_x4& operator-=(v1i_x4& A, v1i_x4 B)
{
    A = A - B;
    return A;
}

inline v1_x4& operator-=(v1_x4& A, v1_x4 B)
{
    A = A - B;
    return A;
}

inline v1u_x4 operator-(v1u_x4 A)
{
    v1u_x4 Result = V1UX4(0) - A;
    return Result;
}

inline v1i_x4 operator-(v1i_x4 A)
{
    v1i_x4 Result = V1IX4(0) - A;
    return Result;
}

inline v1_x4 operator-(v1_x4 A)
{
    v1_x4 Result = V1X4(0) - A;    
    return Result;
}

//
// NOTE: Mul
// 

/*
  // TODO: Int multiplies are not like regular muls, play around here
inline v1i_x4 operator*(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_mul_epi32(A.x, B.x);
#elif MATH_SIMD_ARM
    Result.x = vmulq_s32(A.x, B.x);
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] * B.e[0];
    Result.e[1] = A.e[1] * B.e[1];
    Result.e[2] = A.e[2] * B.e[2];
    Result.e[3] = A.e[3] * B.e[3];
#endif

    return Result;
}
*/

inline v1_x4 operator*(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};
    
#if MATH_SIMD_X64
    Result.x = _mm_mul_ps(A.x, B.x);
#elif MATH_SIMD_ARM
    Result.x = vmulq_f32(A.x, B.x);
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] * B.e[0];
    Result.e[1] = A.e[1] * B.e[1];
    Result.e[2] = A.e[2] * B.e[2];
    Result.e[3] = A.e[3] * B.e[3];
#endif
    
    return Result;
}

/*
inline v1i_x4& operator*=(v1i_x4& A, v1i_x4 B)
{
    A = A * B;
    return A;
}
*/

inline v1_x4& operator*=(v1_x4& A, v1_x4 B)
{
    A = A * B;
    return A;
}

//
// NOTE: Div
// 

inline v1_x4 operator/(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};
    
#if MATH_SIMD_X64
    Result.x = _mm_div_ps(A.x, B.x);
#elif MATH_SIMD_ARM
    Result.x = vdivq_f32(A.x, B.x);
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] / B.e[0];
    Result.e[1] = A.e[1] / B.e[1];
    Result.e[2] = A.e[2] / B.e[2];
    Result.e[3] = A.e[3] / B.e[3];
#endif
    
    return Result;
}

inline v1_x4& operator/=(v1_x4& A, v1_x4 B)
{
    A = A / B;
    return A;
}

//
// NOTE: Logical Not
//

inline v1u_x4 operator~(v1u_x4 V)
{
    v1u_x4 Result = {};
#if MATH_SIMD_X64
    Result.x = _mm_andnot_si128(V.x, _mm_set1_epi32(0xFFFFFFFF));
#elif MATH_SIMD_ARM
    Result.x = vmvnq_u32(V.x);
#elif MATH_SIMD_TEST
    Result.e[0] = ~V.e[0];
    Result.e[1] = ~V.e[1];
    Result.e[2] = ~V.e[2];
    Result.e[3] = ~V.e[3];
#endif
    return Result;
}

inline v1i_x4 operator~(v1i_x4 V)
{
    v1i_x4 Result = {};
#if MATH_SIMD_X64
    Result.x = _mm_andnot_si128(V.x, _mm_set1_epi32(0xFFFFFFFF));
#elif MATH_SIMD_ARM
    Result.x = vmvnq_u32(V.x);
#elif MATH_SIMD_TEST
    Result.e[0] = ~V.e[0];
    Result.e[1] = ~V.e[1];
    Result.e[2] = ~V.e[2];
    Result.e[3] = ~V.e[3];
#endif
    return Result;
}

//
// NOTE: Logical Or
//

inline v1u_x4 operator|(v1u_x4 A, v1u_x4 B)
{
    v1u_x4 Result = {};
    
#if MATH_SIMD_X64
    Result.x = _mm_or_si128(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] | B.e[0];
    Result.e[1] = A.e[1] | B.e[1];
    Result.e[2] = A.e[2] | B.e[2];
    Result.e[3] = A.e[3] | B.e[3];    
#endif

    return Result;
}

inline v1i_x4 operator|(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};
    
#if MATH_SIMD_X64
    Result.x = _mm_or_si128(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] | B.e[0];
    Result.e[1] = A.e[1] | B.e[1];
    Result.e[2] = A.e[2] | B.e[2];
    Result.e[3] = A.e[3] | B.e[3];    
#endif

    return Result;
}

inline v1u_x4& operator|=(v1u_x4& A, v1u_x4 B)
{
    A = A | B;
    return A;
}

inline v1i_x4& operator|=(v1i_x4& A, v1i_x4 B)
{
    A = A | B;
    return A;
}

//
// NOTE: Logical And
//

inline v1u_x4 operator&(v1u_x4 A, v1u_x4 B)
{
    v1u_x4 Result = {};
    
#if MATH_SIMD_X64
    Result.x = _mm_and_si128(A.x, B.x);
#elif MATH_SIMD_ARM
    Result.x = vandq_s32(A.x, B.x);
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] & B.e[0];
    Result.e[1] = A.e[1] & B.e[1];
    Result.e[2] = A.e[2] & B.e[2];
    Result.e[3] = A.e[3] & B.e[3];
#endif

    return Result;
}

inline v1i_x4 operator&(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};
    
#if MATH_SIMD_X64
    Result.x = _mm_and_si128(A.x, B.x);
#elif MATH_SIMD_ARM
    Result.x = vandq_s32(A.x, B.x);
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] & B.e[0];
    Result.e[1] = A.e[1] & B.e[1];
    Result.e[2] = A.e[2] & B.e[2];
    Result.e[3] = A.e[3] & B.e[3];
#endif

    return Result;
}

inline v1_x4 operator&(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};
#if MATH_SIMD_X64
    Result.x = _mm_and_ps(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = ReinterpretF32(ReinterpretU32(A.e[0]) & ReinterpretU32(B.e[0]));
    Result.e[1] = ReinterpretF32(ReinterpretU32(A.e[1]) & ReinterpretU32(B.e[1]));
    Result.e[2] = ReinterpretF32(ReinterpretU32(A.e[2]) & ReinterpretU32(B.e[2]));
    Result.e[3] = ReinterpretF32(ReinterpretU32(A.e[3]) & ReinterpretU32(B.e[3]));
#endif

    return Result;
}

inline v1u_x4& operator&=(v1u_x4& A, v1u_x4 B)
{
    A = A & B;
    return A;
}

inline v1i_x4& operator&=(v1i_x4& A, v1i_x4 B)
{
    A = A & B;
    return A;
}

inline v1_x4& operator&=(v1_x4& A, v1_x4 B)
{
    A = A & B;
    return A;
}

//
// NOTE: Compare Equal
//

inline v1u_x4 operator==(v1u_x4 A, v1u_x4 B)
{
    v1u_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_cmpeq_epi32(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] == B.e[0];
    Result.e[1] = A.e[1] == B.e[1];
    Result.e[2] = A.e[2] == B.e[2];
    Result.e[3] = A.e[3] == B.e[3];
#endif

    return Result;
}

inline v1i_x4 operator==(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_cmpeq_epi32(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] == B.e[0];
    Result.e[1] = A.e[1] == B.e[1];
    Result.e[2] = A.e[2] == B.e[2];
    Result.e[3] = A.e[3] == B.e[3];
#endif

    return Result;
}

inline v1_x4 operator==(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_cmpeq_ps(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] == B.e[0];
    Result.e[1] = A.e[1] == B.e[1];
    Result.e[2] = A.e[2] == B.e[2];
    Result.e[3] = A.e[3] == B.e[3];
#endif

    return Result;
}

//
// NOTE: Compare Not Equal
//

inline v1u_x4 operator!=(v1u_x4 A, v1u_x4 B)
{
    v1u_x4 Result = {};

#if MATH_SIMD_X64
    Result = ~(A == B);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] != B.e[0];
    Result.e[1] = A.e[1] != B.e[1];
    Result.e[2] = A.e[2] != B.e[2];
    Result.e[3] = A.e[3] != B.e[3];
#endif

    return Result;
}

inline v1i_x4 operator!=(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};

#if MATH_SIMD_X64
    Result = ~(A == B);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] != B.e[0];
    Result.e[1] = A.e[1] != B.e[1];
    Result.e[2] = A.e[2] != B.e[2];
    Result.e[3] = A.e[3] != B.e[3];
#endif

    return Result;
}

inline v1_x4 operator!=(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_cmpneq_ps(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] != B.e[0];
    Result.e[1] = A.e[1] != B.e[1];
    Result.e[2] = A.e[2] != B.e[2];
    Result.e[3] = A.e[3] != B.e[3];
#endif

    return Result;
}

//
// NOTE: Compare Less Than
// 

inline v1u_x4 operator<(v1u_x4 A, v1u_x4 B)
{
    v1u_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_cmplt_epi32(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] < B.e[0];
    Result.e[1] = A.e[1] < B.e[1];
    Result.e[2] = A.e[2] < B.e[2];
    Result.e[3] = A.e[3] < B.e[3];
#endif

    return Result;
}

inline v1i_x4 operator<(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_cmplt_epi32(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] < B.e[0];
    Result.e[1] = A.e[1] < B.e[1];
    Result.e[2] = A.e[2] < B.e[2];
    Result.e[3] = A.e[3] < B.e[3];
#endif

    return Result;
}

inline v1_x4 operator<(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_cmplt_ps(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] < B.e[0];
    Result.e[1] = A.e[1] < B.e[1];
    Result.e[2] = A.e[2] < B.e[2];
    Result.e[3] = A.e[3] < B.e[3];
#endif

    return Result;
}

//
// NOTE: Compare Greater Than
// 

inline v1u_x4 operator>(v1u_x4 A, v1u_x4 B)
{
    v1u_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_cmpgt_epi32(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] > B.e[0];
    Result.e[1] = A.e[1] > B.e[1];
    Result.e[2] = A.e[2] > B.e[2];
    Result.e[3] = A.e[3] > B.e[3];
#endif

    return Result;
}

inline v1i_x4 operator>(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_cmpgt_epi32(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] > B.e[0];
    Result.e[1] = A.e[1] > B.e[1];
    Result.e[2] = A.e[2] > B.e[2];
    Result.e[3] = A.e[3] > B.e[3];
#endif

    return Result;
}

inline v1_x4 operator>(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_cmpgt_ps(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] > B.e[0];
    Result.e[1] = A.e[1] > B.e[1];
    Result.e[2] = A.e[2] > B.e[2];
    Result.e[3] = A.e[3] > B.e[3];
#endif

    return Result;
}

//
// NOTE: Compare Less Than or Equal
// 

inline v1u_x4 operator<=(v1u_x4 A, v1u_x4 B)
{
    v1u_x4 Result = {};

#if MATH_SIMD_X64
    Result = ~(A > B);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] <= B.e[0];
    Result.e[1] = A.e[1] <= B.e[1];
    Result.e[2] = A.e[2] <= B.e[2];
    Result.e[3] = A.e[3] <= B.e[3];
#endif

    return Result;
}

inline v1i_x4 operator<=(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};

#if MATH_SIMD_X64
    Result = ~(A > B);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] <= B.e[0];
    Result.e[1] = A.e[1] <= B.e[1];
    Result.e[2] = A.e[2] <= B.e[2];
    Result.e[3] = A.e[3] <= B.e[3];
#endif

    return Result;
}

inline v1_x4 operator<=(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_cmple_ps(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] <= B.e[0];
    Result.e[1] = A.e[1] <= B.e[1];
    Result.e[2] = A.e[2] <= B.e[2];
    Result.e[3] = A.e[3] <= B.e[3];
#endif

    return Result;
}

//
// NOTE: Compare Greater Than or Equal
// 

inline v1u_x4 operator>=(v1u_x4 A, v1u_x4 B)
{
    v1u_x4 Result = {};

#if MATH_SIMD_X64
    Result = ~(A < B);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] >= B.e[0];
    Result.e[1] = A.e[1] >= B.e[1];
    Result.e[2] = A.e[2] >= B.e[2];
    Result.e[3] = A.e[3] >= B.e[3];
#endif

    return Result;
}

inline v1i_x4 operator>=(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};

#if MATH_SIMD_X64
    Result = ~(A < B);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] >= B.e[0];
    Result.e[1] = A.e[1] >= B.e[1];
    Result.e[2] = A.e[2] >= B.e[2];
    Result.e[3] = A.e[3] >= B.e[3];
#endif

    return Result;
}

inline v1_x4 operator>=(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_cmpge_ps(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = A.e[0] >= B.e[0];
    Result.e[1] = A.e[1] >= B.e[1];
    Result.e[2] = A.e[2] >= B.e[2];
    Result.e[3] = A.e[3] >= B.e[3];
#endif

    return Result;
}

// =======================================================================================================================================
// NOTE: Helpful Operators (interface with scalar code cleanly)
// =======================================================================================================================================

//
// NOTE: Add
//

inline v1u_x4 operator+(v1u_x4 A, u32 B)
{
    v1u_x4 Result = A + V1UX4(B);
    return Result;
}

inline v1i_x4 operator+(v1i_x4 A, i32 B)
{
    v1i_x4 Result = A + V1IX4(B);
    return Result;
}

inline v1_x4 operator+(v1_x4 A, f32 B)
{
    v1_x4 Result = A + V1X4(B);
    
    return Result;
}

inline v1u_x4 operator+(u32 A, v1u_x4 B)
{
    v1u_x4 Result = V1UX4(A) + B;
    return Result;
}

inline v1i_x4 operator+(i32 A, v1i_x4 B)
{
    v1i_x4 Result = V1IX4(A) + B;
    return Result;
}

inline v1_x4 operator+(f32 A, v1_x4 B)
{
    v1_x4 Result = V1X4(A) + B;
    return Result;
}

inline v1u_x4& operator+=(v1u_x4& A, u32 B)
{
    A = A + B;
    return A;
}

inline v1i_x4& operator+=(v1i_x4& A, i32 B)
{
    A = A + B;
    return A;
}

inline v1_x4& operator+=(v1_x4& A, f32 B)
{
    A = A + B;
    return A;
}

//
// NOTE: Sub
//

inline v1u_x4 operator-(v1u_x4 A, u32 B)
{
    v1u_x4 Result = A - V1UX4(B);
    return Result;
}

inline v1i_x4 operator-(v1i_x4 A, i32 B)
{
    v1i_x4 Result = A - V1IX4(B);
    return Result;
}

inline v1_x4 operator-(v1_x4 A, f32 B)
{
    v1_x4 Result = A - V1X4(B);
    return Result;
}

inline v1u_x4 operator-(u32 A, v1u_x4 B)
{
    v1u_x4 Result = V1UX4(A) - B;
    return Result;
}

inline v1i_x4 operator-(i32 A, v1i_x4 B)
{
    v1i_x4 Result = V1IX4(A) - B;
    return Result;
}

inline v1_x4 operator-(f32 A, v1_x4 B)
{
    v1_x4 Result = V1X4(A) - B;
    return Result;
}

inline v1u_x4& operator-=(v1u_x4& A, u32 B)
{
    A = A - B;
    return A;
}

inline v1i_x4& operator-=(v1i_x4& A, i32 B)
{
    A = A - B;
    return A;
}

inline v1_x4& operator-=(v1_x4& A, f32 B)
{
    A = A - B;
    return A;
}

//
// NOTE: Mul
//

/*
inline v1i_x4 operator*(v1i_x4 A, i32 B)
{
    v1i_x4 Result = A * V1IX4(B);
    return Result;
}
*/

inline v1_x4 operator*(v1_x4 A, f32 B)
{
    v1_x4 Result = A * V1X4(B);    
    return Result;
}

/*
inline v1i_x4 operator*(i32 A, v1i_x4 B)
{
    v1i_x4 Result = V1IX4(A) * B;
    return Result;
}
*/

inline v1_x4 operator*(f32 A, v1_x4 B)
{
    v1_x4 Result = V1X4(A) * B;    
    return Result;
}

inline v1_x4& operator*=(v1_x4& A, f32 B)
{
    A = A * B;
    return A;
}

/*
inline v1i_x4& operator*=(v1i_x4& A, i32 B)
{
    A = A * B;
    return A;
}
*/

//
// NOTE: Div
//

inline v1_x4 operator/(v1_x4 A, f32 B)
{
    v1_x4 Result = A / V1X4(B);    
    return Result;
}

inline v1_x4 operator/(f32 A, v1_x4 B)
{
    v1_x4 Result = V1X4(A) / B;    
    return Result;
}

inline v1_x4& operator/=(v1_x4& A, f32 B)
{
    A = A / B;
    return A;
}

//
// NOTE: Logical Or
//

inline v1u_x4 operator|(v1u_x4 A, u32 B)
{
    v1u_x4 Result = A | V1UX4(B);    
    return Result;
}

inline v1i_x4 operator|(v1i_x4 A, i32 B)
{
    v1i_x4 Result = A | V1IX4(B);    
    return Result;
}

inline v1u_x4 operator|(u32 A, v1u_x4 B)
{
    v1u_x4 Result = V1UX4(A) | B;    
    return Result;
}

inline v1i_x4 operator|(i32 A, v1i_x4 B)
{
    v1i_x4 Result = V1IX4(A) | B;    
    return Result;
}

inline v1u_x4& operator|=(v1u_x4& A, u32 B)
{
    A = A | B;
    return A;
}

inline v1i_x4& operator|=(v1i_x4& A, i32 B)
{
    A = A | B;
    return A;
}

//
// NOTE: Logical And
//

inline v1u_x4 operator&(v1u_x4 A, u32 B)
{
    v1u_x4 Result = A & V1UX4(B);
    return Result;
}

inline v1i_x4 operator&(v1i_x4 A, i32 B)
{
    v1i_x4 Result = A & V1IX4(B);
    return Result;
}

inline v1_x4 operator&(v1_x4 A, f32 B)
{
    v1_x4 Result = A & V1X4(B);    
    return Result;
}

inline v1u_x4 operator&(u32 A, v1u_x4 B)
{
    v1u_x4 Result = V1UX4(A) & B;
    return Result;
}

inline v1i_x4 operator&(i32 A, v1i_x4 B)
{
    v1i_x4 Result = V1IX4(A) & B;
    return Result;
}

inline v1_x4 operator&(f32 A, v1_x4 B)
{
    v1_x4 Result = V1X4(A) & B;    
    return Result;
}

inline v1u_x4& operator&=(v1u_x4& A, u32 B)
{
    A = A & B;
    return A;
}

inline v1i_x4& operator&=(v1i_x4& A, i32 B)
{
    A = A & B;
    return A;
}

inline v1_x4& operator&=(v1_x4& A, f32 B)
{
    A = A & B;
    return A;
}

// =======================================================================================================================================
// NOTE: Common Math Functions
// =======================================================================================================================================

//
// NOTE: Move Mask
//

inline u32 MoveMask(v1u_x4 V)
{
    // NOTE: Outputs 1 bit for each element if its most sig bit is true
    u32 Result = 0;

#if MATH_SIMD_X64
    Result = _mm_movemask_epi8(V.x);
    Result = ((Result >> 3) & 0x1) | ((Result >> 6) & 0x1) | ((Result >> 10) & 0x1) | ((Result >> 12) & 0x1);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    InvalidCodePath;
#endif

    return Result;
}

inline u32 MoveMask(v1i_x4 V)
{
    // NOTE: Outputs 1 bit for each element if its most sig bit is true
    u32 Result = 0;

#if MATH_SIMD_X64
    Result = _mm_movemask_epi8(V.x);
    Result = ((Result >> 3) & 0x1) | ((Result >> 6) & 0x1) | ((Result >> 10) & 0x1) | ((Result >> 12) & 0x1);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    InvalidCodePath;
#endif

    return Result;
}

inline u32 MoveMask(v1_x4 V)
{
    // NOTE: Outputs 1 bit for each element if its most sig bit is true
    u32 Result = 0;

#if MATH_SIMD_X64
    Result = _mm_movemask_ps(V.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    InvalidCodePath;
#endif

    return Result;
}

//
// NOTE: Min
//

inline v1u_x4 Min(v1u_x4 A, v1u_x4 B)
{
    v1u_x4 Result = {};
    
#if MATH_SIMD_X64
    Result.x = _mm_min_epu32(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = Min(A.e[0], B.e[0]);
    Result.e[1] = Min(A.e[1], B.e[1]);
    Result.e[2] = Min(A.e[2], B.e[2]);
    Result.e[3] = Min(A.e[3], B.e[3]);
#endif
    
    return Result;
}

inline v1i_x4 Min(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};
    
#if MATH_SIMD_X64
    Result.x = _mm_min_epi32(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = Min(A.e[0], B.e[0]);
    Result.e[1] = Min(A.e[1], B.e[1]);
    Result.e[2] = Min(A.e[2], B.e[2]);
    Result.e[3] = Min(A.e[3], B.e[3]);
#endif
    
    return Result;
}

inline v1_x4 Min(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};
    
#if MATH_SIMD_X64
    Result.x = _mm_min_ps(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = Min(A.e[0], B.e[0]);
    Result.e[1] = Min(A.e[1], B.e[1]);
    Result.e[2] = Min(A.e[2], B.e[2]);
    Result.e[3] = Min(A.e[3], B.e[3]);
#endif
    
    return Result;
}

//
// NOTE: Max
//

inline v1u_x4 Max(v1u_x4 A, v1u_x4 B)
{
    v1u_x4 Result = {};
    
#if MATH_SIMD_X64
    Result.x = _mm_max_epu32(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = Max(A.e[0], B.e[0]);
    Result.e[1] = Max(A.e[1], B.e[1]);
    Result.e[2] = Max(A.e[2], B.e[2]);
    Result.e[3] = Max(A.e[3], B.e[3]);
#endif
    
    return Result;
}

inline v1i_x4 Max(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};
    
#if MATH_SIMD_X64
    Result.x = _mm_max_epi32(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = Max(A.e[0], B.e[0]);
    Result.e[1] = Max(A.e[1], B.e[1]);
    Result.e[2] = Max(A.e[2], B.e[2]);
    Result.e[3] = Max(A.e[3], B.e[3]);
#endif
    
    return Result;
}

inline v1_x4 Max(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};
    
#if MATH_SIMD_X64
    Result.x = _mm_max_ps(A.x, B.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = Max(A.e[0], B.e[0]);
    Result.e[1] = Max(A.e[1], B.e[1]);
    Result.e[2] = Max(A.e[2], B.e[2]);
    Result.e[3] = Max(A.e[3], B.e[3]);
#endif
    
    return Result;
}

//
// NOTE: Floor
//

inline v1u_x4 FloorV1UX4(v1_x4 V)
{
    v1u_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_cvtps_epi32((V - 0.5f).x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    InvalidCodePath;
#endif
    
    return Result;
}

inline v1i_x4 FloorV1IX4(v1_x4 V)
{
    v1i_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_cvtps_epi32(V.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    InvalidCodePath;
#endif
    
    return Result;
}

// TODO: This is SSE4.1
inline v1_x4 Floor(v1_x4 A)
{
    v1_x4 Result = {};

#if MATH_SIMD_X64
    Result.x = _mm_floor_ps(A.x);
#elif MATH_SIMD_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    Result.e[0] = FloorF32(A.e[0]);
    Result.e[1] = FloorF32(A.e[1]);
    Result.e[2] = FloorF32(A.e[2]);
    Result.e[3] = FloorF32(A.e[3]);
#endif
    
    return Result;
}


//
// NOTE: Clamp
//

inline v1u_x4 Clamp(v1u_x4 Val, v1u_x4 MinVal, v1u_x4 MaxVal)
{
    v1u_x4 Result = Min(MaxVal, Max(MinVal, Val));
    return Result;
}

inline v1i_x4 Clamp(v1i_x4 Val, v1i_x4 MinVal, v1i_x4 MaxVal)
{
    v1i_x4 Result = Min(MaxVal, Max(MinVal, Val));
    return Result;
}

inline v1_x4 Clamp(v1_x4 Val, v1_x4 MinVal, v1_x4 MaxVal)
{
    v1_x4 Result = Min(MaxVal, Max(MinVal, Val));
    return Result;
}

//
// NOTE: Square Root
//

inline v1_x4 SquareRoot(v1_x4 A)
{
    v1_x4 Result = {};
#if MATH_SIMD_X64
    Result.x = _mm_sqrt_ps(A.x);
#elif MATH_SIMD_ARM
    Result.e[0] = SquareRoot(A.e[0]);
    Result.e[1] = SquareRoot(A.e[1]);
    Result.e[2] = SquareRoot(A.e[2]);
    Result.e[3] = SquareRoot(A.e[3]);
#endif
    return Result;
}

//
// NOTE: Horizontal Min
//

inline u32 HorizontalMin(v1u_x4 V)
{
    u32 Result = 0;

#if MATH_SIMD_X64
    v1u_x4 Shuffle1, Shuffle2;
    Shuffle1.x = _mm_shuffle_epi32(V.x, (0x1 << 0) | (0x0 << 1) | (0x3 << 2) | (0x2 << 3));
    V = Min(V, Shuffle1);
    Shuffle2.x = _mm_shuffle_epi32(V.x, (0x3 << 0) | (0x3 << 1) | (0x0 << 2) | (0x0 << 3));
    V = Min(V, Shuffle2);
    Result = V.e[0];
#elif MATH_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    InvalidCodePath;
#endif

    return Result;
}

inline i32 HorizontalMin(v1i_x4 V)
{
    i32 Result = 0;

#if MATH_SIMD_X64
    v1i_x4 Shuffle1, Shuffle2;
    Shuffle1.x = _mm_shuffle_epi32(V.x, (0x1 << 0) | (0x0 << 1) | (0x3 << 2) | (0x2 << 3));
    V = Min(V, Shuffle1);
    Shuffle2.x = _mm_shuffle_epi32(V.x, (0x3 << 0) | (0x3 << 1) | (0x0 << 2) | (0x0 << 3));
    V = Min(V, Shuffle2);
    Result = V.e[0];
#elif MATH_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    InvalidCodePath;
#endif

    return Result;
}

//
// NOTE: Horizontal Max
//

inline u32 HorizontalMax(v1u_x4 V)
{
    u32 Result = 0;

#if MATH_SIMD_X64
    v1u_x4 Shuffle1, Shuffle2;
    Shuffle1.x = _mm_shuffle_epi32(V.x, (0x1 << 0) | (0x0 << 1) | (0x3 << 2) | (0x2 << 3));
    V = Max(V, Shuffle1);
    Shuffle2.x = _mm_shuffle_epi32(V.x, (0x3 << 0) | (0x3 << 1) | (0x0 << 2) | (0x0 << 3));
    V = Max(V, Shuffle2);
    Result = V.e[0];
#elif MATH_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    InvalidCodePath;
#endif

    return Result;
}

inline i32 HorizontalMax(v1i_x4 V)
{
    i32 Result = 0;

#if MATH_SIMD_X64
    v1i_x4 Shuffle1, Shuffle2;
    Shuffle1.x = _mm_shuffle_epi32(V.x, (0x1 << 0) | (0x0 << 1) | (0x3 << 2) | (0x2 << 3));
    V = Max(V, Shuffle1);
    Shuffle2.x = _mm_shuffle_epi32(V.x, (0x3 << 0) | (0x3 << 1) | (0x0 << 2) | (0x0 << 3));
    V = Max(V, Shuffle2);
    Result = V.e[0];
#elif MATH_ARM
    InvalidCodePath;
#elif MATH_SIMD_TEST
    InvalidCodePath;
#endif

    return Result;
}
