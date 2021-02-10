/*
    NOTE: Every struct has the following possible arguments for its operators:

        - scalar (f32, this gets promoted to the vector type normally)
        - scalar vector (this is like v2, v3, etc, think of it as having a bundle where all 4 vectors are this vector)
        - vector (this is expected behavior)
   
 */

// TODO: Add all comparison operators
// TODO: Add operator[] ops

// =======================================================================================================================================
// NOTE: SSE Helpers
// =======================================================================================================================================

#if MATH_X64

#define SSE_SHUFFLE_MASK(a, b, c, d) ((d << 6) | (c << 4) | (b << 2) | (a << 0))

/*
  NOTE: This is the layout we have

  a00_0, a00_1, a00_2, a00_3,
  a01_0, a01_1, a01_2, a01_3,
  a10_0, a10_1, a10_2, a10_3,
  a11_0, a11_1, a11_2, a11_3,

  This is what we transform into
  
  a00_0, a01_0, a10_0, a11_0, 
  a00_1, a01_1, a10_1, a11_1, 
  a00_2, a01_2, a10_2, a11_2, 
  a00_3, a01_3, a10_3, a11_3, 
         
*/
#define SSE_TRANSPOSE(In1, In2, In3, In4, Out1, Out2, Out3, Out4)       \
    {                                                                   \
        __m128 temp0 = _mm_shuffle_ps(In1, In2, SSE_SHUFFLE_MASK(0, 1, 0, 1)); \
        __m128 temp1 = _mm_shuffle_ps(In1, In2, SSE_SHUFFLE_MASK(2, 3, 2, 3)); \
        __m128 temp2 = _mm_shuffle_ps(In3, In4, SSE_SHUFFLE_MASK(0, 1, 0, 1)); \
        __m128 temp3 = _mm_shuffle_ps(In3, In4, SSE_SHUFFLE_MASK(2, 3, 2, 3)); \
        Out1 = _mm_shuffle_ps(temp0, temp2, SSE_SHUFFLE_MASK(0, 2, 0, 2)); \
        Out2 = _mm_shuffle_ps(temp0, temp2, SSE_SHUFFLE_MASK(1, 3, 1, 3)); \
        Out3 = _mm_shuffle_ps(temp1, temp3, SSE_SHUFFLE_MASK(0, 2, 0, 2)); \
        Out4 = _mm_shuffle_ps(temp1, temp3, SSE_SHUFFLE_MASK(1, 3, 1, 3)); \
    }

#endif

// =======================================================================================================================================
// NOTE: v1u_x4 Simd
// =======================================================================================================================================

// TODO: A lot of this is really just a signed integer, so we have some implicit conversions happening here

//
// NOTE: Init
//

inline v1u_x4 V1UX4(u32 A)
{
    v1u_x4 Result = {};
#if MATH_X64
    Result.x = _mm_set1_epi32(A);
#elif MATH_ARM

    Result.e[0] = A;
    Result.e[1] = A;
    Result.e[2] = A;
    Result.e[3] = A;
    
    //Result.x = vld1q_dup_u32(&A);
#endif
    return Result;
}

inline v1u_x4 V1UX4(void* X)
{
    v1u_x4 Result = {};
#if MATH_X64
    Result.x = _mm_load_si128((__m128i*)X);
#elif MATH_ARM
    //Result.x = vld1q_u32((u32*)X);

    u32* CurrPtr = (u32*)X;
    Result.e[0] = CurrPtr[0];
    Result.e[1] = CurrPtr[1];
    Result.e[2] = CurrPtr[2];
    Result.e[3] = CurrPtr[3];
    
#endif
    return Result;
}

//
// NOTE: Scalar Writes
//

inline void WriteScalar(v1u_x4* V, u32 Index, u32 Value)
{
    Assert(Index < 4);
    V->e[Index] = Value;
}

//
// NOTE: Store
//

inline void StoreAligned(v1u_x4 V, void* Dest)
{
#if MATH_X64
    _mm_store_si128((__m128i*)Dest, V.x);
#elif MATH_ARM
    //vst1q_u32((u32*)Dest, V.x);

    u32* CurrPtr = (u32*)Dest;
    CurrPtr[0] = V.e[0];
    CurrPtr[1] = V.e[1];
    CurrPtr[2] = V.e[2];
    CurrPtr[3] = V.e[3];
    
#endif
}

inline void StoreScalarAligned(v1u_x4 V, void* Dest, mm Stride, u32 NumStored = 4)
{
    // NOTE: We don't have scatter ops so we have to convert registers here
    u8* CurrPtr = (u8*)Dest;

    for (u32 ElementId = 0; ElementId < NumStored; ++ElementId)
    {
        u32* WritePtr = (u32*)CurrPtr;
        WritePtr[0] = V.e[ElementId];
        CurrPtr += Stride;
    }
}

//
// NOTE: Not
//

inline v1u_x4 operator~(v1u_x4 V)
{
    v1u_x4 Result = {};
#if MATH_X64
    Result.x = _mm_andnot_si128(V.x, _mm_set1_epi32(0xFFFFFFFF));
#elif MATH_ARM
    //Result.x = vmvnq_u32(V.x);
    Result.e[0] = ~V.e[0];
    Result.e[1] = ~V.e[1];
    Result.e[2] = ~V.e[2];
    Result.e[3] = ~V.e[3];
#endif
    return Result;
}

//
// NOTE: Comparison
//

inline v1u_x4 operator!=(v1u_x4 V, u32 B)
{
    v1u_x4 Result = {};
#if MATH_X64
    // NOTE: SSE doesn't have a not equal to for ints
    Result.x = _mm_andnot_si128(_mm_cmpeq_epi32(V.x, _mm_set1_epi32(B)), _mm_set1_epi32(0xFFFFFFFF));
#elif MATH_ARM
    // NOTE: ARM doesn't have a not equal to for ints
    //Result.x = vmvnq_u32(vceqq_u32(V.x, vld1q_dup_u32(&B)));
    Result.e[0] = V.e[0] != B;
    Result.e[1] = V.e[1] != B;
    Result.e[2] = V.e[2] != B;
    Result.e[3] = V.e[3] != B;
#endif
    return Result;
}

//
// NOTE: Array indexing
//

inline u32 v1u_x4::operator[](u32 Index)
{
    // IMPORTANT: This is only to be used in places where we are not doing SIMD ops, this is mainly for when
    // we read from memory a vector, and we want to read it as scalars instead of simd registers
    Assert(Index < 4);
    u32 Result = e[Index];
    return Result;
}

// =======================================================================================================================================
// NOTE: v1i_x4 Simd
// =======================================================================================================================================

//
// NOTE: Init
//

inline v1i_x4 V1IX4(i32* X)
{
    v1i_x4 Result = {};
#if MATH_X64
    Result.x = _mm_load_si128((__m128i*)X);
#elif MATH_ARM
    //Result.x = vld1q_s32((i32*)X);
    Result.e[0] = X[0];
    Result.e[1] = X[1];
    Result.e[2] = X[2];
    Result.e[3] = X[3];
#endif
    return Result;
}

inline v1i_x4 V1IX4(i32 X)
{
    v1i_x4 Result = {};
#if MATH_X64
    Result.x = _mm_set1_epi32(X);
#elif MATH_ARM
    //Result.x = vld1q_dup_s32(&X);
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
#if MATH_X64
    Result.x = _mm_set_epi32(W, Z, Y, X);
#elif MATH_ARM
    // TODO: There is probably a better way to do this, but I can't find the intrinsic
    //i32 Temp[4] = { X, Y, Z, W };
    //Result = V1IX4(Temp);

    Result.e[0] = X;
    Result.e[1] = Y;
    Result.e[2] = Z;
    Result.e[3] = W;
    
#endif
    return Result;
}

//
// NOTE: Scalar Writes
//

inline void WriteScalar(v1i_x4* V, u32 Index, i32 Value)
{
    Assert(Index < 4);
    V->e[Index] = Value;
}

//
// NOTE: Store
//

inline void StoreAligned(v1i_x4 V, i32* Dest)
{
#if MATH_X64
    _mm_store_si128((__m128i*)Dest, V.x);
#elif MATH_ARM
    //vst1q_s32((i32*)Dest, V.x);
    Dest[0] = V.e[0];
    Dest[1] = V.e[1];
    Dest[2] = V.e[2];
    Dest[3] = V.e[3];    
#endif
}

inline void StoreScalarAligned(v1i_x4 V, void* Dest, mm Stride, u32 NumStored = 4)
{
    // NOTE: We don't have scatter ops so we have to convert registers here
    u8* CurrPtr = (u8*)Dest;

    for (u32 ElementId = 0; ElementId < NumStored; ++ElementId)
    {
        i32* WritePtr = (i32*)CurrPtr;
        WritePtr[0] = V.e[ElementId];
        CurrPtr += Stride;
    }
}

//
// NOTE: Conversions
//

inline v1i_x4 V1IX4Convert(v1_x4 V)
{
    v1i_x4 Result = {};
#if MATH_X64
    Result.x = _mm_cvtps_epi32(V.x);
#elif MATH_ARM
    //Result.x = vcvtq_s32_f32(V.x);
    Result.e[0] = i32(V.e[0]);
    Result.e[1] = i32(V.e[1]);
    Result.e[2] = i32(V.e[2]);
    Result.e[3] = i32(V.e[3]);    
#endif
    return Result;
}

inline v1i_x4 V1IX4Cast(v1_x4 V)
{
    v1i_x4 Result = {};
#if MATH_X64
    Result.x = _mm_castps_si128(V.x);
#elif MATH_ARM
    //Result.x = vreinterpretq_s32_f32(V.x);
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

inline v1i_x4 operator+(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};
#if MATH_X64
    Result.x = _mm_add_epi32(A.x, B.x);
#elif MATH_ARM
    //Result.x = vaddq_s32(A.x, B.x);
    Result.e[0] = A.e[0] + B.e[0];
    Result.e[1] = A.e[1] + B.e[1];
    Result.e[2] = A.e[2] + B.e[2];
    Result.e[3] = A.e[3] + B.e[3];
#endif

    return Result;
}

inline v1i_x4 operator+(v1i_x4 A, i32 B)
{
    v1i_x4 Result = A + V1IX4(B);
    return Result;
}

inline v1i_x4 operator+(i32 A, v1i_x4 B)
{
    v1i_x4 Result = V1IX4(A) + B;
    return Result;
}

inline v1i_x4& operator+=(v1i_x4& A, i32 B)
{
    A = A + B;
    return A;
}

inline v1i_x4& operator+=(v1i_x4& A, v1i_x4 B)
{
    A = A + B;
    return A;
}

//
// NOTE: Sub
//

inline v1i_x4 operator-(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};
#if MATH_X64
    Result.x = _mm_sub_epi32(A.x, B.x);
#elif MATH_ARM
    //Result.x = vsubq_s32(A.x, B.x);
    Result.e[0] = A.e[0] - B.e[0];
    Result.e[1] = A.e[1] - B.e[1];
    Result.e[2] = A.e[2] - B.e[2];
    Result.e[3] = A.e[3] - B.e[3];
#endif
    return Result;
}

inline v1i_x4 operator-(v1i_x4 A, i32 B)
{
    v1i_x4 Result = A - V1IX4(B);
    return Result;
}

inline v1i_x4 operator-(i32 A, v1i_x4 B)
{
    v1i_x4 Result = V1IX4(A) - B;
    return Result;
}

inline v1i_x4& operator-=(v1i_x4& A, i32 B)
{
    A = A - B;
    return A;
}

inline v1i_x4& operator-=(v1i_x4& A, v1i_x4 B)
{
    A = A - B;
    return A;
}

//
// NOTE: Negation
//

inline v1i_x4 operator-(v1i_x4 A)
{
    v1i_x4 Result = 0 - A;
    return Result;
}

//
// NOTE: Mul
//

inline v1i_x4 operator*(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};
#if MATH_X64
    Result.x = _mm_mul_epi32(A.x, B.x);
#elif MATH_ARM
    //Result.x = vmulq_s32(A.x, B.x);
    Result.e[0] = A.e[0] * B.e[0];
    Result.e[1] = A.e[1] * B.e[1];
    Result.e[2] = A.e[2] * B.e[2];
    Result.e[3] = A.e[3] * B.e[3];
#endif

    return Result;
}

inline v1i_x4 operator*(v1i_x4 A, i32 B)
{
    v1i_x4 Result = A * V1IX4(B);
    return Result;
}

inline v1i_x4 operator*(i32 A, v1i_x4 B)
{
    v1i_x4 Result = V1IX4(A) * B;
    return Result;
}

inline v1i_x4& operator*=(v1i_x4& A, v1i_x4 B)
{
    A = A * B;
    return A;
}

inline v1i_x4& operator*=(v1i_x4& A, i32 B)
{
    A = A * B;
    return A;
}

//
// NOTE: Not
//

inline v1i_x4 operator~(v1i_x4 A)
{
    v1i_x4 Result = {};
#if MATH_X64
    Result.x = _mm_andnot_si128(A.x, _mm_set1_epi32(0xFFFFFFFF));
#elif MATH_ARM
    Result.e[0] = ~A.e[0];
    Result.e[1] = ~A.e[1];
    Result.e[2] = ~A.e[2];
    Result.e[3] = ~A.e[3];
#endif
    return Result;
}

//
// NOTE: And
//

inline v1i_x4 operator&(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};
#if MATH_X64
    Result.x = _mm_and_si128(A.x, B.x);
#elif MATH_ARM
    //Result.x = vandq_s32(A.x, B.x);
    Result.e[0] = A.e[0] & B.e[0];
    Result.e[1] = A.e[1] & B.e[1];
    Result.e[2] = A.e[2] & B.e[2];
    Result.e[3] = A.e[3] & B.e[3];
#endif

    return Result;
}

inline v1i_x4 operator&(v1i_x4 A, i32 B)
{
    v1i_x4 Result = A & V1IX4(B);
    return Result;
}

inline v1i_x4 operator&(i32 A, v1i_x4 B)
{
    v1i_x4 Result = V1IX4(A) & B;
    return Result;
}

inline v1i_x4& operator&=(v1i_x4& A, v1i_x4 B)
{
    A = A & B;
    return A;
}

inline v1i_x4& operator&=(v1i_x4& A, i32 B)
{
    A = A & B;
    return A;
}

//
// NOTE: Or
//

inline v1i_x4 operator|(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};
#if MATH_X64
    Result.x = _mm_or_si128(A.x, B.x);
#elif MATH_ARM
    Result.e[0] = A.e[0] | B.e[0];
    Result.e[1] = A.e[1] | B.e[1];
    Result.e[2] = A.e[2] | B.e[2];
    Result.e[3] = A.e[3] | B.e[3];    
#endif

    return Result;
}

inline v1i_x4 operator|(v1i_x4 A, i32 B)
{
    v1i_x4 Result = A | V1IX4(B);    
    return Result;
}

inline v1i_x4 operator|(i32 A, v1i_x4 B)
{
    v1i_x4 Result = V1IX4(A) | B;    
    return Result;
}

inline v1i_x4& operator|=(v1i_x4& A, v1i_x4 B)
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
// NOTE: Xor
//

inline v1i_x4 operator^(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};
#if MATH_X64
    Result.x = _mm_xor_si128(A.x, B.x);
#elif MATH_ARM
    Result.e[0] = A.e[0] ^ B.e[0];
    Result.e[1] = A.e[1] ^ B.e[1];
    Result.e[2] = A.e[2] ^ B.e[2];
    Result.e[3] = A.e[3] ^ B.e[3];
#endif

    return Result;
}

inline v1i_x4 operator^(v1i_x4 A, i32 B)
{
    v1i_x4 Result = A ^ V1IX4(B);
    return Result;
}

inline v1i_x4 operator^(i32 A, v1i_x4 B)
{
    v1i_x4 Result = V1IX4(A) ^ B;
    return Result;
}

inline v1i_x4& operator^=(v1i_x4& A, i32 B)
{
    A = A ^ B;
    return A;
}

inline v1i_x4& operator^=(v1i_x4& A, v1i_x4 B)
{
    A = A ^ B;
    return A;
}

//
// NOTE: Shift Left
//

inline v1i_x4 operator<<(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};
#if MATH_X64
    Result.x = _mm_sll_epi32(A.x, B.x);
#elif MATH_ARM
    //Result.x = vshlq_s32(A.x, B.x);
    Result.e[0] = A.e[0] << B.e[0];
    Result.e[1] = A.e[1] << B.e[1];
    Result.e[2] = A.e[2] << B.e[2];
    Result.e[3] = A.e[3] << B.e[3];
#endif

    return Result;
}

inline v1i_x4 operator<<(v1i_x4 A, i32 B)
{
    v1i_x4 Result = {};
#if MATH_X64
    Result.x = _mm_slli_epi32(A.x, B);
#elif MATH_ARM
    // TODO: There is a shift by integer intrinsic but it needs to be a constant int
    //Result.x = vshlq_s32(A.x, vld1q_dup_s32(&B));
    Result = A << V1IX4(B);
#endif

    return Result;
}

inline v1i_x4& operator<<=(v1i_x4& A, i32 B)
{
    A = A << B;
    return A;
}

inline v1i_x4& operator<<=(v1i_x4& A, v1i_x4 B)
{
    A = A << B;
    return A;
}

//
// NOTE: Shift Right
//

inline v1i_x4 operator>>(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};
#if MATH_X64
    Result.x = _mm_srl_epi32(A.x, B.x);
#elif MATH_ARM
    // TODO: Shift right intrinsics seem to be limited, idk why
    //Result.x = vshrq_s32(A.x, B.x);
    Result.e[0] = A.e[0] >> B.e[0];
    Result.e[1] = A.e[1] >> B.e[1];
    Result.e[2] = A.e[2] >> B.e[2];
    Result.e[3] = A.e[3] >> B.e[3];
#endif

    return Result;
}

inline v1i_x4 operator>>(v1i_x4 A, i32 B)
{
    v1i_x4 Result = {};
#if MATH_X64
    Result.x = _mm_srli_epi32(A.x, B);
#elif MATH_ARM
    // TODO: Shift right intrinsics seem to be limited, idk why
    Result = A >> V1IX4(B);
#endif

    return Result;
}

inline v1i_x4& operator>>=(v1i_x4& A, i32 B)
{
    A = A >> B;
    return A;
}

inline v1i_x4& operator>>=(v1i_x4& A, v1i_x4 B)
{
    A = A >> B;
    return A;
}

//
// NOTE: Compare Equal
//

inline v1i_x4 operator==(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};
#if MATH_X64
    Result.x = _mm_cmpeq_epi32(A.x, B.x);
#elif MATH_ARM
    //Result.x = vceqq_s32(A.x, B.x);
    Result.e[0] = A.e[0] == B.e[0];
    Result.e[1] = A.e[1] == B.e[1];
    Result.e[2] = A.e[2] == B.e[2];
    Result.e[3] = A.e[3] == B.e[3];
#endif

    return Result;
}

inline v1i_x4 operator==(v1i_x4 A, i32 B)
{
    v1i_x4 Result = A == V1IX4(B);
    return Result;
}

inline v1i_x4 operator==(i32 A, v1i_x4 B)
{
    v1i_x4 Result = V1IX4(A) == B;
    return Result;
}

// =======================================================================================================================================
// NOTE: v1_x4 Simd
// =======================================================================================================================================

//
// NOTE: Init
//

inline v1_x4 V1X4(f32* X)
{
    v1_x4 Result = {};
#if MATH_X64
    Result.x = _mm_load_ps(X);
#elif MATH_ARM
    //Result.x = vld1q_f32(X);
    Result.e[0] = X[0];
    Result.e[1] = X[1];
    Result.e[2] = X[2];
    Result.e[3] = X[3];
#endif
    return Result;
}

inline v1_x4 V1X4(f32 X)
{
    v1_x4 Result = {};
#if MATH_X64
    Result.x = _mm_set1_ps(X);
#elif MATH_ARM
    //Result.x = vld1q_dup_f32(&X);
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
#if MATH_X64
    Result.x = _mm_set_ps(W, Z, Y, X);
#elif MATH_ARM
    // TODO: Is there a intrinsic for this?
    Result.e[0] = X;
    Result.e[1] = Y;
    Result.e[2] = Z;
    Result.e[3] = W;
#endif
    return Result;
}

//
// NOTE: Scalar Writes
//

inline void WriteScalar(v1_x4* V, u32 Index, f32 Value)
{
    Assert(Index < 4);
    V->e[Index] = Value;
}

//
// NOTE: Store
//

inline void StoreAligned(v1_x4 V, f32* Dest)
{
#if MATH_X64
    _mm_store_ps(Dest, V.x);
#elif MATH_ARM
    //vst1q_f32((f32*)Dest, V.x);
    Dest[0] = V.e[0];
    Dest[1] = V.e[1];
    Dest[2] = V.e[2];
    Dest[3] = V.e[3];
#endif
}

//
// NOTE: Conversions
//

inline v1_x4 V1X4Convert(v1i_x4 V)
{
    v1_x4 Result = {};
#if MATH_X64
    Result.x = _mm_cvtepi32_ps(V.x);
#elif MATH_ARM
    // TODO: Is there a intrinsic for this?
    Result.e[0] = f32(V.e[0]);
    Result.e[1] = f32(V.e[1]);
    Result.e[2] = f32(V.e[2]);
    Result.e[3] = f32(V.e[3]);
#endif
    return Result;
}

inline v1_x4 V1X4Cast(v1i_x4 V)
{
    v1_x4 Result = {};
#if MATH_X64
    Result.x = _mm_castsi128_ps(V.x);
#elif MATH_ARM
    //Result.x = vreinterpretq_f32_s32(V.x);
    Result.e[0] = ReinterpretF32(V.e[0]);
    Result.e[1] = ReinterpretF32(V.e[1]);
    Result.e[2] = ReinterpretF32(V.e[2]);
    Result.e[3] = ReinterpretF32(V.e[3]);
#endif
    return Result;
}

inline v1_x4 V1X4Cast(v1u_x4 V)
{
    v1_x4 Result = {};
#if MATH_X64
    Result.x = _mm_castsi128_ps(V.x);
#elif MATH_ARM
    //Result.x = vreinterpretq_f32_u32(V.x);
    Result.e[0] = ReinterpretF32(V.e[0]);
    Result.e[1] = ReinterpretF32(V.e[1]);
    Result.e[2] = ReinterpretF32(V.e[2]);
    Result.e[3] = ReinterpretF32(V.e[3]);
#endif
    return Result;
}

//
// NOTE: Add
//

inline v1_x4 operator+(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};
#if MATH_X64
    Result.x = _mm_add_ps(A.x, B.x);
#elif MATH_ARM
    //Result.x = vaddq_f32(A.x, B.x);
    Result.e[0] = A.e[0] + B.e[0];
    Result.e[1] = A.e[1] + B.e[1];
    Result.e[2] = A.e[2] + B.e[2];
    Result.e[3] = A.e[3] + B.e[3];
#endif
    
    return Result;
}

inline v1_x4 operator+(v1_x4 A, f32 B)
{
    v1_x4 Result = A + V1X4(B);
    
    return Result;
}

inline v1_x4 operator+(f32 A, v1_x4 B)
{
    v1_x4 Result = V1X4(A) + B;
    return Result;
}

inline v1_x4& operator+=(v1_x4& A, f32 B)
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

inline v1_x4 operator-(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};
#if MATH_X64
    Result.x = _mm_sub_ps(A.x, B.x);
#elif MATH_ARM
    //Result.x = vsubq_f32(A.x, B.x);
    Result.e[0] = A.e[0] - B.e[0];
    Result.e[1] = A.e[1] - B.e[1];
    Result.e[2] = A.e[2] - B.e[2];
    Result.e[3] = A.e[3] - B.e[3];
#endif
    
    return Result;
}

inline v1_x4 operator-(v1_x4 A, f32 B)
{
    v1_x4 Result = A - V1X4(B);
    return Result;
}

inline v1_x4 operator-(f32 A, v1_x4 B)
{
    v1_x4 Result = V1X4(A) - B;
    return Result;
}

inline v1_x4& operator-=(v1_x4& A, f32 B)
{
    A = A - B;
    return A;
}

inline v1_x4& operator-=(v1_x4& A, v1_x4 B)
{
    A = A - B;
    return A;
}

//
// NOTE: Negation
//

inline v1_x4 operator-(v1_x4 A)
{
    v1_x4 Result = 0 - A;    
    return Result;
}

//
// NOTE: Mul
//

inline v1_x4 operator*(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};
#if MATH_X64
    Result.x = _mm_mul_ps(A.x, B.x);
#elif MATH_ARM
    //Result.x = vmulq_f32(A.x, B.x);
    Result.e[0] = A.e[0] * B.e[0];
    Result.e[1] = A.e[1] * B.e[1];
    Result.e[2] = A.e[2] * B.e[2];
    Result.e[3] = A.e[3] * B.e[3];
#endif
    
    return Result;
}

inline v1_x4 operator*(v1_x4 A, f32 B)
{
    v1_x4 Result = A * V1X4(B);    
    return Result;
}

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
#if MATH_X64
    Result.x = _mm_div_ps(A.x, B.x);
#elif MATH_ARM
    //Result.x = vdivq_f32(A.x, B.x);
    Result.e[0] = A.e[0] / B.e[0];
    Result.e[1] = A.e[1] / B.e[1];
    Result.e[2] = A.e[2] / B.e[2];
    Result.e[3] = A.e[3] / B.e[3];
#endif
    
    return Result;
}

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

inline v1_x4& operator/=(v1_x4& A, v1_x4 B)
{
    A = A / B;
    return A;
}

//
// NOTE: Not
//

inline v1_x4 operator~(v1_x4 A)
{
    v1_x4 Result = {};
#if MATH_X64
    Result.x = _mm_andnot_ps(A.x, _mm_set1_ps(0xFFFFFFFF));
#elif MATH_ARM
    Result.e[0] = ReinterpretF32(~ReinterpretU32(A.e[0]));
    Result.e[1] = ReinterpretF32(~ReinterpretU32(A.e[1]));
    Result.e[2] = ReinterpretF32(~ReinterpretU32(A.e[2]));
    Result.e[3] = ReinterpretF32(~ReinterpretU32(A.e[3]));
#endif
    return Result;
}

//
// NOTE: And
//

inline v1_x4 operator&(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};
#if MATH_X64
    Result.x = _mm_and_ps(A.x, B.x);
#elif MATH_ARM
    Result.e[0] = ReinterpretF32(ReinterpretU32(A.e[0]) & ReinterpretU32(B.e[0]));
    Result.e[1] = ReinterpretF32(ReinterpretU32(A.e[1]) & ReinterpretU32(B.e[1]));
    Result.e[2] = ReinterpretF32(ReinterpretU32(A.e[2]) & ReinterpretU32(B.e[2]));
    Result.e[3] = ReinterpretF32(ReinterpretU32(A.e[3]) & ReinterpretU32(B.e[3]));
#endif

    return Result;
}

inline v1_x4 operator&(v1_x4 A, f32 B)
{
    v1_x4 Result = A & V1X4(B);    
    return Result;
}

inline v1_x4 operator&(f32 A, v1_x4 B)
{
    v1_x4 Result = V1X4(A) & B;    
    return Result;
}

inline v1_x4& operator&=(v1_x4& A, f32 B)
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
// NOTE: Or
//

inline v1_x4 operator|(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};
#if MATH_X64
    Result.x = _mm_or_ps(A.x, B.x);
#elif MATH_ARM
    Result.e[0] = ReinterpretF32(ReinterpretU32(A.e[0]) | ReinterpretU32(B.e[0]));
    Result.e[1] = ReinterpretF32(ReinterpretU32(A.e[1]) | ReinterpretU32(B.e[1]));
    Result.e[2] = ReinterpretF32(ReinterpretU32(A.e[2]) | ReinterpretU32(B.e[2]));
    Result.e[3] = ReinterpretF32(ReinterpretU32(A.e[3]) | ReinterpretU32(B.e[3]));
#endif

    return Result;
}

inline v1_x4 operator|(v1_x4 A, f32 B)
{
    v1_x4 Result = A | V1X4(B);    
    return Result;
}

inline v1_x4 operator|(f32 A, v1_x4 B)
{
    v1_x4 Result = V1X4(A) | B;    
    return Result;
}

inline v1_x4& operator|=(v1_x4& A, f32 B)
{
    A = A | B;
    return A;
}

inline v1_x4& operator|=(v1_x4& A, v1_x4 B)
{
    A = A | B;
    return A;
}

//
// NOTE: Xor
//

inline v1_x4 operator^(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};
#if MATH_X64
    Result.x = _mm_xor_ps(A.x, B.x);
#elif MATH_ARM
    Result.e[0] = ReinterpretF32(ReinterpretU32(A.e[0]) ^ ReinterpretU32(B.e[0]));
    Result.e[1] = ReinterpretF32(ReinterpretU32(A.e[1]) ^ ReinterpretU32(B.e[1]));
    Result.e[2] = ReinterpretF32(ReinterpretU32(A.e[2]) ^ ReinterpretU32(B.e[2]));
    Result.e[3] = ReinterpretF32(ReinterpretU32(A.e[3]) ^ ReinterpretU32(B.e[3]));
#endif

    return Result;
}

inline v1_x4 operator^(v1_x4 A, f32 B)
{
    v1_x4 Result = A ^ V1X4(B);    
    return Result;
}

inline v1_x4 operator^(f32 A, v1_x4 B)
{
    v1_x4 Result = V1X4(A) ^ B;    
    return Result;
}

inline v1_x4& operator^=(v1_x4& A, f32 B)
{
    A = A ^ B;
    return A;
}

inline v1_x4& operator^=(v1_x4& A, v1_x4 B)
{
    A = A ^ B;
    return A;
}

//
// NOTE: Mask Ops
//

inline v1_x4 MaskedWrite(v1_x4 A, v1_x4 B, v1u_x4 Mask, v1u_x4 NotMask)
{
    v1_x4 Result = {};
    Result = (A & V1X4Cast(Mask)) | (A & V1X4Cast(NotMask));
    return Result;
}

inline v1_x4 MaskedWrite(v1_x4 A, v1_x4 B, v1u_x4 Mask)
{
    v1_x4 Result = MaskedWrite(A, B, Mask, ~Mask);
    return Result;
}

// =======================================================================================================================================
// NOTE: v2_x4
// =======================================================================================================================================

//
// NOTE: Init
//

inline v2_x4 V2X4(f32* X, f32* Y)
{
    v2_x4 Result = {};
    Result.x = V1X4(X);
    Result.y = V1X4(Y);
    return Result;
}

inline v2_x4 V2X4(f32 X, f32 Y)
{
    v2_x4 Result = {};
    Result.x = V1X4(X);
    Result.y = V1X4(Y);
    return Result;
}

inline v2_x4 V2X4(v2 V)
{
    v2_x4 Result = {};
    Result.x = V1X4(V.x);
    Result.y = V1X4(V.y);
    return Result;
}

inline v2_x4 V2X4(v1_x4 X, v1_x4 Y)
{
    v2_x4 Result = {};
    Result.x = X;
    Result.y = Y;
    return Result;
}

inline v2_x4 V2X4(v2_soa V, u32 Id)
{
    v2_x4 Result = V2X4(V.x + Id, V.y + Id);
    return Result;
}

//
// NOTE: Scalar Writes
//

inline void WriteScalar(v2_x4* V, u32 Index, v2 Value)
{
    WriteScalar(&V->x, Index, Value.x);
    WriteScalar(&V->y, Index, Value.y);
}

//
// NOTE: Store
//

inline void StoreAligned(v2_x4 V, f32* DestX, f32* DestY)
{
    StoreAligned(V.x, DestX);
    StoreAligned(V.y, DestY);
}

inline void StoreAligned(v2_x4 V, v2_soa Dest, u32 Id)
{
    Assert((Id % 4) == 0);
    StoreAligned(V.x, Dest.x + Id);
    StoreAligned(V.y, Dest.y + Id);
}

//
// NOTE: Add
//

inline v2_x4 operator+(v2_x4 A, v2_x4 B)
{
    v2_x4 Result = {};
    Result.x = A.x + B.x;
    Result.y = A.y + B.y;
    
    return Result;
}

inline v2_x4 operator+(v2_x4 A, f32 B)
{
    v2_x4 Result = A + V2X4(B, B);
    return Result;
}

inline v2_x4 operator+(f32 A, v2_x4 B)
{
    v2_x4 Result = V2X4(A, A) + B;
    return Result;
}

inline v2_x4& operator+=(v2_x4& A, f32 B)
{
    A = A + B;
    return A;
}

inline v2_x4& operator+=(v2_x4& A, v2_x4 B)
{
    A = A + B;
    return A;
}

//
// NOTE: Sub
//

inline v2_x4 operator-(v2_x4 A, v2_x4 B)
{
    v2_x4 Result = {};
    Result.x = A.x - B.x;
    Result.y = A.y - B.y;
    
    return Result;
}

inline v2_x4 operator-(v2_x4 A, f32 B)
{
    v2_x4 Result = A - V2X4(B, B);
    return Result;
}

inline v2_x4 operator-(f32 A, v2_x4 B)
{
    v2_x4 Result = V2X4(A, A) - B;
    return Result;
}

inline v2_x4 operator-(v2 A, v2_x4 B)
{
    v2_x4 Result = V2X4(A) - B;
    return Result;
}

inline v2_x4 operator-(v2_x4 A, v2 B)
{
    v2_x4 Result = A - V2X4(B);
    return Result;
}

inline v2_x4& operator-=(v2_x4& A, f32 B)
{
    A = A - B;
    return A;
}

inline v2_x4& operator-=(v2_x4& A, v2 B)
{
    A = A - B;
    return A;
}

inline v2_x4& operator-=(v2_x4& A, v2_x4 B)
{
    A = A - B;
    return A;
}

//
// NOTE: Negation
//

inline v2_x4 operator-(v2_x4 A)
{
    v2_x4 Result = V2X4(0.0f, 0.0f) - A;
    return Result;
}

//
// NOTE: Mul
//

inline v2_x4 operator*(v2_x4 A, f32 B)
{
    v2_x4 Result = {};
    Result.x = A.x * B;
    Result.y = A.y * B;
    
    return Result;
}

inline v2_x4 operator*(f32 A, v2_x4 B)
{
    v2_x4 Result = {};
    Result.x = A * B.x;
    Result.y = A * B.y;
    
    return Result;
}

inline v2_x4& operator*=(v2_x4& A, f32 B)
{
    A = A * B;
    return A;
}

inline v2_x4 operator*(v2_x4 A, v2 B)
{
    v2_x4 Result = {};
    Result.x = A.x * B.x;
    Result.y = A.y * B.y;
    
    return Result;
}

inline v2_x4 operator*(v2 A, v2_x4 B)
{
    v2_x4 Result = {};
    Result.x = A.x * B.x;
    Result.y = A.y * B.y;
    
    return Result;
}

inline v2_x4& operator*=(v2_x4& A, v2 B)
{
    A = A * B;
    return A;
}

inline v2_x4 operator*(v2_x4 A, v1_x4 B)
{
    v2_x4 Result = {};
    Result.x = A.x * B;
    Result.y = A.y * B;
    
    return Result;
}

inline v2_x4& operator*=(v2_x4& A, v1_x4 B)
{
    A = A * B;
    return A;
}

inline v2_x4 operator*(v1_x4 B, v2_x4 A)
{
    v2_x4 Result = {};
    Result.x = A.x * B;
    Result.y = A.y * B;
    
    return Result;
}

inline v2_x4 operator*(v2_x4 A, v2_x4 B)
{
    v2_x4 Result = {};
    Result.x = A.x * B.x;
    Result.y = A.y * B.y;
    
    return Result;
}

inline v2_x4& operator*=(v2_x4& A, v2_x4 B)
{
    A = A * B;
    return A;
}

//
// NOTE: Div
//

inline v2_x4 operator/(v2_x4 A, f32 B)
{
    v2_x4 Result = {};
    Result.x = A.x / B;
    Result.y = A.y / B;
    
    return Result;
}

inline v2_x4 operator/(f32 A, v2_x4 B)
{
    v2_x4 Result = {};
    Result.x = A / B.x;
    Result.y = A / B.y;
    
    return Result;
}

inline v2_x4& operator/=(v2_x4& A, f32 B)
{
    A = A / B;
    return A;
}

inline v2_x4 operator/(v2_x4 A, v2 B)
{
    v2_x4 Result = {};
    Result.x = A.x / B.x;
    Result.y = A.y / B.y;
    
    return Result;
}

inline v2_x4 operator/(v2 A, v2_x4 B)
{
    v2_x4 Result = {};
    Result.x = A.x / B.x;
    Result.y = A.y / B.y;
    
    return Result;
}

inline v2_x4& operator/=(v2_x4& A, v2 B)
{
    A = A / B;
    return A;
}

inline v2_x4 operator/(v2_x4 A, v1_x4 B)
{
    v2_x4 Result = {};
    Result.x = A.x / B;
    Result.y = A.y / B;
    
    return Result;
}

inline v2_x4& operator/=(v2_x4& A, v1_x4 B)
{
    A = A / B;
    return A;
}

inline v2_x4 operator/(v2_x4 A, v2_x4 B)
{
    v2_x4 Result = {};
    Result.x = A.x / B.x;
    Result.y = A.y / B.y;
    
    return Result;
}

inline v2_x4& operator/=(v2_x4& A, v2_x4 B)
{
    A = A / B;
    return A;
}

//
// NOTE: And
//

inline v2_x4 operator&(v2_x4 A, f32 B)
{
    v2_x4 Result = {};
    Result.x = A.x & B;
    Result.y = A.y & B;
    
    return Result;
}

inline v2_x4 operator&(f32 A, v2_x4 B)
{
    v2_x4 Result = {};
    Result.x = A & B.x;
    Result.y = A & B.y;
    
    return Result;
}

inline v2_x4& operator&=(v2_x4& A, f32 B)
{
    A = A & B;
    return A;
}

inline v2_x4 operator&(v2_x4 A, v2_x4 B)
{
    v2_x4 Result = {};
    Result.x = A.x & B.x;
    Result.y = A.y & B.y;
    return Result;
}

inline v2_x4& operator&=(v2_x4& A, v2_x4 B)
{
    A = A & B;
    return A;
}

//
// NOTE: Or
//

inline v2_x4 operator|(v2_x4 A, f32 B)
{
    v2_x4 Result = {};
    Result.x = A.x | B;
    Result.y = A.y | B;
    
    return Result;
}

inline v2_x4 operator|(f32 A, v2_x4 B)
{
    v2_x4 Result = {};
    Result.x = A | B.x;
    Result.y = A | B.y;
    
    return Result;
}

inline v2_x4& operator|=(v2_x4& A, f32 B)
{
    A = A | B;
    return A;
}

inline v2_x4 operator|(v2_x4 A, v2_x4 B)
{
    v2_x4 Result = {};
    Result.x = A.x | B.x;
    Result.y = A.y | B.y;
    return Result;
}

inline v2_x4& operator|=(v2_x4& A, v2_x4 B)
{
    A = A | B;
    return A;
}

//
// NOTE: Xor
//

inline v2_x4 operator^(v2_x4 A, f32 B)
{
    v2_x4 Result = {};
    Result.x = A.x ^ B;
    Result.y = A.y ^ B;
    
    return Result;
}

inline v2_x4 operator^(f32 A, v2_x4 B)
{
    v2_x4 Result = {};
    Result.x = A ^ B.x;
    Result.y = A ^ B.y;
    
    return Result;
}

inline v2_x4& operator^=(v2_x4& A, f32 B)
{
    A = A ^ B;
    return A;
}

inline v2_x4 operator^(v2_x4 A, v2_x4 B)
{
    v2_x4 Result = {};
    Result.x = A.x ^ B.x;
    Result.y = A.y ^ B.y;
    return Result;
}

inline v2_x4& operator^=(v2_x4& A, v2_x4 B)
{
    A = A ^ B;
    return A;
}

//
// NOTE: Mask Ops
//

inline v2_x4 MaskedWrite(v2_x4 A, v2_x4 B, v1u_x4 Mask, v1u_x4 NotMask)
{
    v2_x4 Result = {};
    Result.x = MaskedWrite(A.x, B.x, Mask, NotMask);
    Result.y = MaskedWrite(A.y, B.y, Mask, NotMask);
    return Result;
}

inline v2_x4 MaskedWrite(v2_x4 A, v2_x4 B, v1u_x4 Mask)
{
    v2_x4 Result = MaskedWrite(A, B, Mask, ~Mask);
    return Result;
}

// =======================================================================================================================================
// NOTE: v3_x4
// =======================================================================================================================================

//
// NOTE: Init
//

inline v3_x4 V3X4(f32* X, f32* Y, f32* Z)
{
    v3_x4 Result = {};
    Result.x = V1X4(X);
    Result.y = V1X4(Y);
    Result.z = V1X4(Z);
    return Result;
}

inline v3_x4 V3X4(f32 X, f32 Y, f32 Z)
{
    v3_x4 Result = {};
    Result.x = V1X4(X);
    Result.y = V1X4(Y);
    Result.z = V1X4(Z);
    return Result;
}

inline v3_x4 V3X4(v3 V)
{
    v3_x4 Result = {};
    Result.x = V1X4(V.x);
    Result.y = V1X4(V.y);
    Result.z = V1X4(V.z);
    return Result;
}

inline v3_x4 V3X4(v1_x4 X, v1_x4 Y, v1_x4 Z)
{
    v3_x4 Result = {};
    Result.x = X;
    Result.y = Y;
    Result.z = Z;
    return Result;
}

inline v3_x4 V3X4(v2_x4 XY, v1_x4 Z)
{
    v3_x4 Result = {};
    Result.x = XY.x;
    Result.y = XY.y;
    Result.z = Z;
    return Result;
}

inline v3_x4 V3X4(v3_soa V, u32 Id)
{
    v3_x4 Result = V3X4(V.x + Id, V.y + Id, V.z + Id);
    return Result;
}

//
// NOTE: Scalar Writes
//

inline void WriteScalar(v3_x4* V, u32 Index, v3 Value)
{
    WriteScalar(&V->x, Index, Value.x);
    WriteScalar(&V->y, Index, Value.y);
    WriteScalar(&V->z, Index, Value.z);
}

//
// NOTE: Store
//

inline void StoreAligned(v3_x4 V, f32* DestX, f32* DestY, f32* DestZ)
{
    StoreAligned(V.x, DestX);
    StoreAligned(V.y, DestY);
    StoreAligned(V.z, DestZ);
}

inline void StoreAligned(v3_x4 V, v3_soa Dest, u32 Id)
{
    Assert((Id % 4) == 0);
    StoreAligned(V.x, Dest.x + Id);
    StoreAligned(V.y, Dest.y + Id);
    StoreAligned(V.z, Dest.z + Id);
}

//
// NOTE: Add
//

inline v3_x4 operator+(v3_x4 A, f32 B)
{
    v3_x4 Result = {};
    Result.x = A.x + B;
    Result.y = A.y + B;
    Result.z = A.z + B;
    
    return Result;
}

inline v3_x4 operator+(f32 A, v3_x4 B)
{
    v3_x4 Result = {};
    Result.x = A + B.x;
    Result.y = A + B.y;
    Result.z = A + B.z;
    
    return Result;
}

inline v3_x4& operator+=(v3_x4& A, f32 B)
{
    A = A + B;
    return A;
}

inline v3_x4 operator+(v3_x4 A, v3_x4 B)
{
    v3_x4 Result = {};
    Result.x = A.x + B.x;
    Result.y = A.y + B.y;
    Result.z = A.z + B.z;
    
    return Result;
}

inline v3_x4& operator+=(v3_x4& A, v3_x4 B)
{
    A = A + B;
    return A;
}

//
// NOTE: Negation
//

inline v3_x4 operator-(v3_x4 A)
{
    v3_x4 Result = {};
    Result.x = V1X4(0.0f) - A.x;
    Result.y = V1X4(0.0f) - A.y;
    Result.z = V1X4(0.0f) - A.z;
    
    return Result;
}

//
// NOTE: Sub
//

inline v3_x4 operator-(v3_x4 A, f32 B)
{
    v3_x4 Result = {};
    Result.x = A.x - B;
    Result.y = A.y - B;
    Result.z = A.z - B;
    
    return Result;
}

inline v3_x4 operator-(f32 A, v3_x4 B)
{
    v3_x4 Result = {};
    Result.x = A - B.x;
    Result.y = A - B.y;
    Result.z = A - B.z;
    
    return Result;
}

inline v3_x4& operator-=(v3_x4& A, f32 B)
{
    A = A - B;
    return A;
}

inline v3_x4 operator-(v3_x4 A, v3 B)
{
    v3_x4 Result = {};
    Result.x = A.x - B.x;
    Result.y = A.y - B.y;
    Result.z = A.z - B.z;
    
    return Result;
}

inline v3_x4 operator-(v3 A, v3_x4 B)
{
    v3_x4 Result = {};
    Result.x = A.x - B.x;
    Result.y = A.y - B.y;
    Result.z = A.z - B.z;
    
    return Result;
}

inline v3_x4& operator-=(v3_x4& A, v3 B)
{
    A = A - B;
    return A;
}

inline v3_x4 operator-(v3_x4 A, v3_x4 B)
{
    v3_x4 Result = {};
    Result.x = A.x - B.x;
    Result.y = A.y - B.y;
    Result.z = A.z - B.z;
    
    return Result;
}

inline v3_x4& operator-=(v3_x4& A, v3_x4 B)
{
    A = A - B;
    return A;
}

//
// NOTE: Mul
//

inline v3_x4 operator*(v3_x4 A, f32 B)
{
    v3_x4 Result = {};
    Result.x = A.x * B;
    Result.y = A.y * B;
    Result.z = A.z * B;
    
    return Result;
}

inline v3_x4 operator*(f32 A, v3_x4 B)
{
    v3_x4 Result = {};
    Result.x = A * B.x;
    Result.y = A * B.y;
    Result.z = A * B.z;
    
    return Result;
}

inline v3_x4& operator*=(v3_x4& A, f32 B)
{
    A = A * B;
    return A;
}

inline v3_x4 operator*(v3_x4 A, v3 B)
{
    v3_x4 Result = {};
    Result.x = A.x * B.x;
    Result.y = A.y * B.y;
    Result.z = A.z * B.z;
    
    return Result;
}

inline v3_x4 operator*(v3 A, v3_x4 B)
{
    v3_x4 Result = {};
    Result.x = A.x * B.x;
    Result.y = A.y * B.y;
    Result.z = A.z * B.z;
    
    return Result;
}

inline v3_x4& operator*=(v3_x4& A, v3 B)
{
    A = A * B;
    return A;
}

inline v3_x4 operator*(v3_x4 A, v1_x4 B)
{
    v3_x4 Result = {};
    Result.x = A.x * B;
    Result.y = A.y * B;
    Result.z = A.z * B;
    
    return Result;
}

inline v3_x4& operator*=(v3_x4& A, v1_x4 B)
{
    A = A * B;
    return A;
}

inline v3_x4 operator*(v1_x4 B, v3_x4 A)
{
    v3_x4 Result = {};
    Result.x = A.x * B;
    Result.y = A.y * B;
    Result.z = A.z * B;
    
    return Result;
}

inline v3_x4 operator*(v3_x4 A, v3_x4 B)
{
    v3_x4 Result = {};
    Result.x = A.x * B.x;
    Result.y = A.y * B.y;
    Result.z = A.z * B.z;
    
    return Result;
}

inline v3_x4& operator*=(v3_x4& A, v3_x4 B)
{
    A = A * B;
    return A;
}

//
// NOTE: Div
//

inline v3_x4 operator/(v3_x4 A, f32 B)
{
    v3_x4 Result = {};
    Result.x = A.x / B;
    Result.y = A.y / B;
    Result.z = A.z / B;
    
    return Result;
}

inline v3_x4 operator/(f32 A, v3_x4 B)
{
    v3_x4 Result = {};
    Result.x = A / B.x;
    Result.y = A / B.y;
    Result.z = A / B.z;
    
    return Result;
}

inline v3_x4& operator/=(v3_x4& A, f32 B)
{
    A = A / B;
    return A;
}

inline v3_x4 operator/(v3_x4 A, v3 B)
{
    v3_x4 Result = {};
    Result.x = A.x / B.x;
    Result.y = A.y / B.y;
    Result.z = A.z / B.z;
    
    return Result;
}

inline v3_x4 operator/(v3 A, v3_x4 B)
{
    v3_x4 Result = {};
    Result.x = A.x / B.x;
    Result.y = A.y / B.y;
    Result.z = A.z / B.z;
    
    return Result;
}

inline v3_x4& operator/=(v3_x4& A, v3 B)
{
    A = A / B;
    return A;
}

inline v3_x4 operator/(v3_x4 A, v1_x4 B)
{
    v3_x4 Result = {};
    Result.x = A.x / B;
    Result.y = A.y / B;
    Result.z = A.z / B;
    
    return Result;
}

inline v3_x4& operator/=(v3_x4& A, v1_x4 B)
{
    A = A / B;
    return A;
}

inline v3_x4 operator/(v3_x4 A, v3_x4 B)
{
    v3_x4 Result = {};
    Result.x = A.x / B.x;
    Result.y = A.y / B.y;
    Result.z = A.z / B.z;
    
    return Result;
}

inline v3_x4& operator/=(v3_x4& A, v3_x4 B)
{
    A = A / B;
    return A;
}

//
// NOTE: And
//

inline v3_x4 operator&(v3_x4 A, f32 B)
{
    v3_x4 Result = {};
    Result.x = A.x & B;
    Result.y = A.y & B;
    Result.z = A.z & B;
    
    return Result;
}

inline v3_x4 operator&(f32 A, v3_x4 B)
{
    v3_x4 Result = {};
    Result.x = A & B.x;
    Result.y = A & B.y;
    Result.z = A & B.z;
    
    return Result;
}

inline v3_x4& operator&=(v3_x4& A, f32 B)
{
    A = A & B;
    return A;
}

inline v3_x4 operator&(v3_x4 A, v3_x4 B)
{
    v3_x4 Result = {};
    Result.x = A.x & B.x;
    Result.y = A.y & B.y;
    Result.z = A.z & B.z;
    return Result;
}

inline v3_x4& operator&=(v3_x4& A, v3_x4 B)
{
    A = A & B;
    return A;
}

//
// NOTE: Or
//

inline v3_x4 operator|(v3_x4 A, f32 B)
{
    v3_x4 Result = {};
    Result.x = A.x | B;
    Result.y = A.y | B;
    Result.z = A.z | B;
    
    return Result;
}

inline v3_x4 operator|(f32 A, v3_x4 B)
{
    v3_x4 Result = {};
    Result.x = A | B.x;
    Result.y = A | B.y;
    Result.z = A | B.z;
    
    return Result;
}

inline v3_x4& operator|=(v3_x4& A, f32 B)
{
    A = A | B;
    return A;
}

inline v3_x4 operator|(v3_x4 A, v3_x4 B)
{
    v3_x4 Result = {};
    Result.x = A.x | B.x;
    Result.y = A.y | B.y;
    Result.z = A.z | B.z;
    return Result;
}

inline v3_x4& operator|=(v3_x4& A, v3_x4 B)
{
    A = A | B;
    return A;
}

//
// NOTE: Xor
//

inline v3_x4 operator^(v3_x4 A, f32 B)
{
    v3_x4 Result = {};
    Result.x = A.x ^ B;
    Result.y = A.y ^ B;
    Result.z = A.z ^ B;
    
    return Result;
}

inline v3_x4 operator^(f32 A, v3_x4 B)
{
    v3_x4 Result = {};
    Result.x = A ^ B.x;
    Result.y = A ^ B.y;
    Result.z = A ^ B.z;
    
    return Result;
}

inline v3_x4& operator^=(v3_x4& A, f32 B)
{
    A = A ^ B;
    return A;
}

inline v3_x4 operator^(v3_x4 A, v3_x4 B)
{
    v3_x4 Result = {};
    Result.x = A.x ^ B.x;
    Result.y = A.y ^ B.y;
    Result.z = A.z ^ B.z;
    return Result;
}

inline v3_x4& operator^=(v3_x4& A, v3_x4 B)
{
    A = A ^ B;
    return A;
}

//
// NOTE: Mask Ops
//

inline v3_x4 MaskedWrite(v3_x4 A, v3_x4 B, v1u_x4 Mask, v1u_x4 NotMask)
{
    v3_x4 Result = {};
    Result.x = MaskedWrite(A.x, B.x, Mask, NotMask);
    Result.y = MaskedWrite(A.y, B.y, Mask, NotMask);
    Result.z = MaskedWrite(A.z, B.z, Mask, NotMask);
    return Result;
}

inline v3_x4 MaskedWrite(v3_x4 A, v3_x4 B, v1u_x4 Mask)
{
    v3_x4 Result = MaskedWrite(A, B, Mask, ~Mask);
    return Result;
}

// =======================================================================================================================================
// NOTE: v4_x4
// =======================================================================================================================================

//
// NOTE: Init
//

inline v4_x4 V4X4(f32* X, f32* Y, f32* Z, f32* W)
{
    v4_x4 Result = {};
    Result.x = V1X4(X);
    Result.y = V1X4(Y);
    Result.z = V1X4(Z);
    Result.w = V1X4(W);
    return Result;
}

inline v4_x4 V4X4(f32 X, f32 Y, f32 Z, f32 W)
{
    v4_x4 Result = {};
    Result.x = V1X4(X);
    Result.y = V1X4(Y);
    Result.z = V1X4(Z);
    Result.w = V1X4(W);
    return Result;
}

inline v4_x4 V4X4(v4 V)
{
    v4_x4 Result = {};
    Result.x = V1X4(V.x);
    Result.y = V1X4(V.y);
    Result.z = V1X4(V.z);
    Result.w = V1X4(V.w);
    return Result;
}

inline v4_x4 V4X4(v1_x4 X, v1_x4 Y, v1_x4 Z, v1_x4 W)
{
    v4_x4 Result = {};
    Result.x = X;
    Result.y = Y;
    Result.z = Z;
    Result.w = W;
    return Result;
}

inline v4_x4 V4X4(v3_x4 XYZ, v1_x4 W)
{
    v4_x4 Result = {};
    Result.x = XYZ.x;
    Result.y = XYZ.y;
    Result.z = XYZ.z;
    Result.w = W;
    return Result;
}

inline v4_x4 V4X4(v4_soa V, u32 Id)
{
    v4_x4 Result = V4X4(V.x + Id, V.y + Id, V.z + Id, V.w + Id);
    return Result;
}

//
// NOTE: Scalar Writes
//

inline void WriteScalar(v4_x4* V, u32 Index, v4 Value)
{
    WriteScalar(&V->x, Index, Value.x);
    WriteScalar(&V->y, Index, Value.y);
    WriteScalar(&V->z, Index, Value.z);
    WriteScalar(&V->w, Index, Value.w);
}

//
// NOTE: Store
//

inline void StoreAligned(v4_x4 V, f32* DestX, f32* DestY, f32* DestZ, f32* DestW)
{
    StoreAligned(V.x, DestX);
    StoreAligned(V.y, DestY);
    StoreAligned(V.z, DestZ);
    StoreAligned(V.w, DestW);
}

inline void StoreAligned(v4_x4 V, v4_soa Dest, u32 Id)
{
    Assert((Id % 4) == 0);
    StoreAligned(V.x, Dest.x + Id);
    StoreAligned(V.y, Dest.y + Id);
    StoreAligned(V.z, Dest.z + Id);
    StoreAligned(V.w, Dest.w + Id);
}

inline void StoreScalarAligned(v4_x4 V, void* Dest, mm Stride, u32 NumStored = 4)
{
#if MATH_X64
    __m128 OutVec[4];
    SSE_TRANSPOSE(V.x.x, V.y.x, V.z.x, V.w.x, OutVec[0], OutVec[1], OutVec[2], OutVec[3]);

    u8* CurrPtr = (u8*)Dest;
    for (u32 ElementId = 0; ElementId < NumStored; ++ElementId)
    {
        _mm_store_ps((f32*)CurrPtr, OutVec[ElementId]);
        CurrPtr += Stride;
    }
#elif MATH_ARM
    u8* CurrPtr = (u8*)Dest;
    for (u32 ElementId = 0; ElementId < NumStored; ++ElementId)
    {
        f32* WriteFloat = (f32*)CurrPtr;
        WriteFloat[0] = V.x.e[ElementId];
        WriteFloat[1] = V.y.e[ElementId];
        WriteFloat[2] = V.z.e[ElementId];
        WriteFloat[3] = V.w.e[ElementId];
        CurrPtr += Stride;
    }
#endif
}

//
// NOTE: Add
//

inline v4_x4 operator+(v4_x4 A, f32 B)
{
    v4_x4 Result = {};
    Result.x = A.x + B;
    Result.y = A.y + B;
    Result.z = A.z + B;
    Result.w = A.w + B;
    
    return Result;
}

inline v4_x4 operator+(f32 A, v4_x4 B)
{
    v4_x4 Result = {};
    Result.x = A + B.x;
    Result.y = A + B.y;
    Result.z = A + B.z;
    Result.w = A + B.w;
    
    return Result;
}

inline v4_x4& operator+=(v4_x4& A, f32 B)
{
    A = A + B;
    return A;
}

inline v4_x4 operator+(v4_x4 A, v4_x4 B)
{
    v4_x4 Result = {};
    Result.x = A.x + B.x;
    Result.y = A.y + B.y;
    Result.z = A.z + B.z;
    Result.w = A.w + B.w;
    
    return Result;
}

inline v4_x4& operator+=(v4_x4& A, v4_x4 B)
{
    A = A + B;
    return A;
}

//
// NOTE: Negation
//

inline v4_x4 operator-(v4_x4 A)
{
    v4_x4 Result = {};
    Result.x = V1X4(0.0f) - A.x;
    Result.y = V1X4(0.0f) - A.y;
    Result.z = V1X4(0.0f) - A.z;
    Result.w = V1X4(0.0f) - A.w;
    
    return Result;
}

//
// NOTE: Sub
//

inline v4_x4 operator-(v4_x4 A, f32 B)
{
    v4_x4 Result = {};
    Result.x = A.x - B;
    Result.y = A.y - B;
    Result.z = A.z - B;
    Result.w = A.w - B;
    
    return Result;
}

inline v4_x4 operator-(f32 A, v4_x4 B)
{
    v4_x4 Result = {};
    Result.x = A - B.x;
    Result.y = A - B.y;
    Result.z = A - B.z;
    Result.w = A - B.w;
    
    return Result;
}

inline v4_x4& operator-=(v4_x4& A, f32 B)
{
    A = A - B;
    return A;
}

inline v4_x4 operator-(v4_x4 A, v4 B)
{
    v4_x4 Result = {};
    Result.x = A.x - B.x;
    Result.y = A.y - B.y;
    Result.z = A.z - B.z;
    Result.w = A.w - B.w;
    
    return Result;
}

inline v4_x4 operator-(v4 A, v4_x4 B)
{
    v4_x4 Result = {};
    Result.x = A.x - B.x;
    Result.y = A.y - B.y;
    Result.z = A.z - B.z;
    Result.w = A.w - B.w;
    
    return Result;
}

inline v4_x4& operator-=(v4_x4& A, v4 B)
{
    A = A - B;
    return A;
}

inline v4_x4 operator-(v4_x4 A, v4_x4 B)
{
    v4_x4 Result = {};
    Result.x = A.x - B.x;
    Result.y = A.y - B.y;
    Result.z = A.z - B.z;
    Result.w = A.w - B.w;
    
    return Result;
}

inline v4_x4& operator-=(v4_x4& A, v4_x4 B)
{
    A = A - B;
    return A;
}

//
// NOTE: Mul
//

inline v4_x4 operator*(v4_x4 A, f32 B)
{
    v4_x4 Result = {};
    Result.x = A.x * B;
    Result.y = A.y * B;
    Result.z = A.z * B;
    Result.w = A.w * B;
    
    return Result;
}

inline v4_x4 operator*(f32 A, v4_x4 B)
{
    v4_x4 Result = {};
    Result.x = A * B.x;
    Result.y = A * B.y;
    Result.z = A * B.z;
    Result.w = A * B.w;
    
    return Result;
}

inline v4_x4& operator*=(v4_x4& A, f32 B)
{
    A = A * B;
    return A;
}

inline v4_x4 operator*(v4_x4 A, v4 B)
{
    v4_x4 Result = {};
    Result.x = A.x * B.x;
    Result.y = A.y * B.y;
    Result.z = A.z * B.z;
    Result.w = A.w * B.w;
    
    return Result;
}

inline v4_x4 operator*(v4 A, v4_x4 B)
{
    v4_x4 Result = {};
    Result.x = A.x * B.x;
    Result.y = A.y * B.y;
    Result.z = A.z * B.z;
    Result.w = A.w * B.w;
    
    return Result;
}

inline v4_x4& operator*=(v4_x4& A, v4 B)
{
    A = A * B;
    return A;
}

inline v4_x4 operator*(v4_x4 A, v1_x4 B)
{
    v4_x4 Result = {};
    Result.x = A.x * B;
    Result.y = A.y * B;
    Result.z = A.z * B;
    Result.w = A.w * B;
    
    return Result;
}

inline v4_x4& operator*=(v4_x4& A, v1_x4 B)
{
    A = A * B;
    return A;
}

inline v4_x4 operator*(v1_x4 B, v4_x4 A)
{
    v4_x4 Result = {};
    Result.x = A.x * B;
    Result.y = A.y * B;
    Result.z = A.z * B;
    Result.w = A.w * B;
    
    return Result;
}

inline v4_x4 operator*(v4_x4 A, v4_x4 B)
{
    v4_x4 Result = {};
    Result.x = A.x * B.x;
    Result.y = A.y * B.y;
    Result.z = A.z * B.z;
    Result.w = A.w * B.w;
    
    return Result;
}

inline v4_x4& operator*=(v4_x4& A, v4_x4 B)
{
    A = A * B;
    return A;
}

//
// NOTE: Div
//

inline v4_x4 operator/(v4_x4 A, f32 B)
{
    v4_x4 Result = {};
    Result.x = A.x / B;
    Result.y = A.y / B;
    Result.z = A.z / B;
    Result.w = A.w / B;
    
    return Result;
}

inline v4_x4 operator/(f32 A, v4_x4 B)
{
    v4_x4 Result = {};
    Result.x = A / B.x;
    Result.y = A / B.y;
    Result.z = A / B.z;
    Result.w = A / B.w;
    
    return Result;
}

inline v4_x4& operator/=(v4_x4& A, f32 B)
{
    A = A / B;
    return A;
}

inline v4_x4 operator/(v4_x4 A, v4 B)
{
    v4_x4 Result = {};
    Result.x = A.x / B.x;
    Result.y = A.y / B.y;
    Result.z = A.z / B.z;
    Result.w = A.w / B.w;
    
    return Result;
}

inline v4_x4 operator/(v4 A, v4_x4 B)
{
    v4_x4 Result = {};
    Result.x = A.x / B.x;
    Result.y = A.y / B.y;
    Result.z = A.z / B.z;
    Result.w = A.w / B.w;
    
    return Result;
}

inline v4_x4& operator/=(v4_x4& A, v4 B)
{
    A = A / B;
    return A;
}

inline v4_x4 operator/(v4_x4 A, v1_x4 B)
{
    v4_x4 Result = {};
    Result.x = A.x / B;
    Result.y = A.y / B;
    Result.z = A.z / B;
    Result.w = A.w / B;
    
    return Result;
}

inline v4_x4& operator/=(v4_x4& A, v1_x4 B)
{
    A = A / B;
    return A;
}

inline v4_x4 operator/(v4_x4 A, v4_x4 B)
{
    v4_x4 Result = {};
    Result.x = A.x / B.x;
    Result.y = A.y / B.y;
    Result.z = A.z / B.z;
    Result.w = A.w / B.w;
    
    return Result;
}

inline v4_x4& operator/=(v4_x4& A, v4_x4 B)
{
    A = A / B;
    return A;
}

//
// NOTE: And
//

inline v4_x4 operator&(v4_x4 A, f32 B)
{
    v4_x4 Result = {};
    Result.x = A.x & B;
    Result.y = A.y & B;
    Result.z = A.z & B;
    Result.w = A.w & B;
    
    return Result;
}

inline v4_x4 operator&(f32 A, v4_x4 B)
{
    v4_x4 Result = {};
    Result.x = A & B.x;
    Result.y = A & B.y;
    Result.z = A & B.z;
    Result.w = A & B.w;
    
    return Result;
}

inline v4_x4& operator&=(v4_x4& A, f32 B)
{
    A = A & B;
    return A;
}

inline v4_x4 operator&(v4_x4 A, v4_x4 B)
{
    v4_x4 Result = {};
    Result.x = A.x & B.x;
    Result.y = A.y & B.y;
    Result.z = A.z & B.z;
    Result.w = A.w & B.w;
    return Result;
}

inline v4_x4& operator&=(v4_x4& A, v4_x4 B)
{
    A = A & B;
    return A;
}

//
// NOTE: Or
//

inline v4_x4 operator|(v4_x4 A, f32 B)
{
    v4_x4 Result = {};
    Result.x = A.x | B;
    Result.y = A.y | B;
    Result.z = A.z | B;
    Result.w = A.w | B;
    
    return Result;
}

inline v4_x4 operator|(f32 A, v4_x4 B)
{
    v4_x4 Result = {};
    Result.x = A | B.x;
    Result.y = A | B.y;
    Result.z = A | B.z;
    Result.w = A | B.w;
    
    return Result;
}

inline v4_x4& operator|=(v4_x4& A, f32 B)
{
    A = A | B;
    return A;
}

inline v4_x4 operator|(v4_x4 A, v4_x4 B)
{
    v4_x4 Result = {};
    Result.x = A.x | B.x;
    Result.y = A.y | B.y;
    Result.z = A.z | B.z;
    Result.w = A.w | B.w;
    return Result;
}

inline v4_x4& operator|=(v4_x4& A, v4_x4 B)
{
    A = A | B;
    return A;
}

//
// NOTE: Xor
//

inline v4_x4 operator^(v4_x4 A, f32 B)
{
    v4_x4 Result = {};
    Result.x = A.x ^ B;
    Result.y = A.y ^ B;
    Result.z = A.z ^ B;
    Result.w = A.w ^ B;
    
    return Result;
}

inline v4_x4 operator^(f32 A, v4_x4 B)
{
    v4_x4 Result = {};
    Result.x = A ^ B.x;
    Result.y = A ^ B.y;
    Result.z = A ^ B.z;
    Result.w = A ^ B.w;
    
    return Result;
}

inline v4_x4& operator^=(v4_x4& A, f32 B)
{
    A = A ^ B;
    return A;
}

inline v4_x4 operator^(v4_x4 A, v4_x4 B)
{
    v4_x4 Result = {};
    Result.x = A.x ^ B.x;
    Result.y = A.y ^ B.y;
    Result.z = A.z ^ B.z;
    Result.w = A.w ^ B.w;
    return Result;
}

inline v4_x4& operator^=(v4_x4& A, v4_x4 B)
{
    A = A ^ B;
    return A;
}

//
// NOTE: Color spaces
//

inline v4_x4 SRGBToLinear(v4_x4 Texel)
{
    v4_x4 Result = (1.0f / 255.0f)*Texel;
    
    Result.r = Square(Result.r);
    Result.g = Square(Result.g);
    Result.b = Square(Result.b);

    return Result;
}

inline v4_x4 LinearToSRGB(v4_x4 Texel)
{
    v4_x4 Result = Texel;
    
    Result.r = SquareRoot(Texel.r);
    Result.g = SquareRoot(Texel.g);
    Result.b = SquareRoot(Texel.b);

    Result *= 255.0f;

    return Result;
}

inline v4_x4 PreMulAlpha(v4_x4 Texel)
{
    v4_x4 Result = Texel;
    Result.xyz *= Texel.a;

    return Result;
}

//
// NOTE: Mask Ops
//

inline v4_x4 MaskedWrite(v4_x4 A, v4_x4 B, v1u_x4 Mask, v1u_x4 NotMask)
{
    v4_x4 Result = {};
    Result.x = MaskedWrite(A.x, B.x, Mask, NotMask);
    Result.y = MaskedWrite(A.y, B.y, Mask, NotMask);
    Result.z = MaskedWrite(A.z, B.z, Mask, NotMask);
    Result.z = MaskedWrite(A.w, B.w, Mask, NotMask);
    return Result;
}

inline v4_x4 MaskedWrite(v4_x4 A, v4_x4 B, v1u_x4 Mask)
{
    v4_x4 Result = MaskedWrite(A, B, Mask, ~Mask);
    return Result;
}

// =======================================================================================================================================
// NOTE: q4_x4
// =======================================================================================================================================

//
// NOTE: Init
//

inline q4_x4 Q4X4(f32* X, f32* Y, f32* Z, f32* W)
{
    q4_x4 Result = {};
    Result.x = V1X4(X);
    Result.y = V1X4(Y);
    Result.z = V1X4(Z);
    Result.w = V1X4(W);
    return Result;
}

inline q4_x4 Q4X4(f32 X, f32 Y, f32 Z, f32 W)
{
    q4_x4 Result = {};
    Result.x = V1X4(X);
    Result.y = V1X4(Y);
    Result.z = V1X4(Z);
    Result.w = V1X4(W);
    return Result;
}

inline q4_x4 Q4X4(q4 V)
{
    q4_x4 Result = {};
    Result.x = V1X4(V.x);
    Result.y = V1X4(V.y);
    Result.z = V1X4(V.z);
    Result.w = V1X4(V.w);
    return Result;
}

inline q4_x4 Q4X4(q4_soa V, u32 Id)
{
    q4_x4 Result = Q4X4(V.x + Id, V.y + Id, V.z + Id, V.w + Id);
    return Result;
}

inline q4_x4 Q4AxisAngle(v3_x4 Axis, v1_x4 Angle)
{
    q4_x4 Result = {};
    Result.xyz = Axis*Sin(V1X4(0.5f)*Angle);
    Result.w = Cos(V1X4(0.5f)*Angle);

    return Result;
}

inline q4_x4 Q4EulerAngles(v1_x4 Yaw, v1_x4 Pitch, v1_x4 Roll)
{
    // NOTE: https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
    v1_x4 CosYaw = Cos(Yaw*V1X4(0.5f));
    v1_x4 SinYaw = Sin(Yaw*V1X4(0.5f));
    v1_x4 CosPitch = Cos(Pitch*V1X4(0.5f));
    v1_x4 SinPitch = Sin(Pitch*V1X4(0.5f));
    v1_x4 CosRoll = Cos(Roll*V1X4(0.5f));
    v1_x4 SinRoll = Sin(Roll*V1X4(0.5f));

    q4_x4 Result = {};
    Result.x = CosRoll * SinYaw * CosPitch - SinRoll * CosYaw * SinPitch;
    Result.y = CosRoll * CosYaw * SinPitch + SinRoll * SinYaw * CosPitch;
    Result.z = SinRoll * CosYaw * CosPitch - CosRoll * SinYaw * SinPitch;
    Result.w = CosRoll * CosYaw * CosPitch + SinRoll * SinYaw * SinPitch;

    return Result;
}

//
// NOTE: Scalar Writes
//

inline void WriteScalar(q4_x4* V, u32 Index, q4 Value)
{
    WriteScalar(&V->x, Index, Value.x);
    WriteScalar(&V->y, Index, Value.y);
    WriteScalar(&V->z, Index, Value.z);
    WriteScalar(&V->w, Index, Value.w);
}

//
// NOTE: Store
//

inline void StoreAligned(q4_x4 Q, f32* DestX, f32* DestY, f32* DestZ, f32* DestW)
{
    StoreAligned(Q.x, DestX);
    StoreAligned(Q.y, DestY);
    StoreAligned(Q.z, DestZ);
    StoreAligned(Q.w, DestW);
}

inline void StoreAligned(q4_x4 Q, q4_soa Dest, u32 Id)
{
    Assert((Id % 4) == 0);
    StoreAligned(Q.x, Dest.x + Id);
    StoreAligned(Q.y, Dest.y + Id);
    StoreAligned(Q.z, Dest.z + Id);
    StoreAligned(Q.w, Dest.w + Id);
}

//
// NOTE: Mul
//

inline q4_x4 operator*(q4_x4 A, v1_x4 B)
{
    q4_x4 Result = {};
    Result.x = A.x * B;
    Result.y = A.y * B;
    Result.z = A.z * B;
    Result.w = A.w * B;
    
    return Result;
}

inline q4_x4& operator*=(q4_x4& A, v1_x4 B)
{
    A = A * B;
    return A;
}

inline q4_x4 operator*(v1_x4 B, q4_x4 A)
{
    q4_x4 Result = {};
    Result.x = A.x * B;
    Result.y = A.y * B;
    Result.z = A.z * B;
    Result.w = A.w * B;
    
    return Result;
}

inline q4_x4 operator*(q4_x4 A, q4_x4 B)
{
    q4_x4 Result = {};
    Result.x = A.x * B.w + A.y * B.z - A.z * B.y + A.w * B.x;
    Result.y = -A.x * B.z + A.y * B.w + A.z * B.x + A.w * B.y;
    Result.z = A.x * B.y - A.y * B.x + A.z * B.w + A.w * B.z;
    Result.w = -A.x * B.x - A.y * B.y - A.z * B.z + A.w * B.w;
    
    return Result;
}

inline q4_x4& operator*=(q4_x4& A, q4_x4 B)
{
    A = A * B;
    return A;
}

//
// NOTE: Div
//

inline q4_x4 operator/(q4_x4 A, v1_x4 B)
{
    q4_x4 Result = {};
    Result.x = A.x / B;
    Result.y = A.y / B;
    Result.z = A.z / B;
    Result.w = A.w / B;
    
    return Result;
}

inline q4_x4& operator/=(q4_x4& A, v1_x4 B)
{
    A = A / B;
    return A;
}

//
// NOTE: Quaternion Functions
//

inline q4_x4 Conjugate(q4_x4 Q)
{
    q4_x4 Result = {};
    Result.x = -Q.x;
    Result.y = -Q.y;
    Result.z = -Q.z;
    Result.w = Q.w;

    return Result;
}

inline q4_x4 Inverse(q4_x4 Q)
{
    q4_x4 Result = Conjugate(Q) / LengthSquared(Q);
    return Result;
}

inline m3_x4 Q4ToM3(q4_x4 Q)
{
    m3_x4 Result = {};
    Result.v[0].x = 1.0f - 2.0f*Square(Q.y) - 2.0f*Square(Q.z);
    Result.v[0].y = 2.0f*Q.x*Q.y + 2.0f*Q.z*Q.w;
    Result.v[0].z = 2.0f*Q.x*Q.z - 2.0f*Q.y*Q.w;

    Result.v[1].x = 2.0f*Q.x*Q.y - 2.0f*Q.z*Q.w;
    Result.v[1].y = 1.0f - 2.0f*Square(Q.x) - 2.0f*Square(Q.z);
    Result.v[1].z = 2.0f*Q.y*Q.z + 2.0f*Q.x*Q.w;

    Result.v[2].x = 2.0f*Q.x*Q.z + 2.0f*Q.y*Q.w;
    Result.v[2].y = 2.0f*Q.y*Q.z - 2.0f*Q.x*Q.w;
    Result.v[2].z = 1.0f - 2.0f*Square(Q.x) - 2.0f*Square(Q.y);

    return Result;
}

inline m4_x4 Q4ToM4(q4_x4 Q)
{
    m4_x4 Result = {};
    Result.v[0].x = 1.0f - 2.0f*Square(Q.y) - 2.0f*Square(Q.z);
    Result.v[0].y = 2.0f*Q.x*Q.y + 2.0f*Q.z*Q.w;
    Result.v[0].z = 2.0f*Q.x*Q.z - 2.0f*Q.y*Q.w;
    Result.v[0].w = V1X4(0.0f);

    Result.v[1].x = 2.0f*Q.x*Q.y - 2.0f*Q.z*Q.w;
    Result.v[1].y = 1.0f - 2.0f*Square(Q.x) - 2.0f*Square(Q.z);
    Result.v[1].z = 2.0f*Q.y*Q.z + 2.0f*Q.x*Q.w;
    Result.v[1].w = V1X4(0.0f);

    Result.v[2].x = 2.0f*Q.x*Q.z + 2.0f*Q.y*Q.w;
    Result.v[2].y = 2.0f*Q.y*Q.z - 2.0f*Q.x*Q.w;
    Result.v[2].z = 1.0f - 2.0f*Square(Q.x) - 2.0f*Square(Q.y);
    Result.v[2].w = V1X4(0.0f);

    Result.v[3].x = V1X4(0.0f);
    Result.v[3].y = V1X4(0.0f);
    Result.v[3].z = V1X4(0.0f);
    Result.v[3].w = V1X4(1.0f);

    return Result;
}

//
// NOTE: Mask Ops
//

inline q4_x4 MaskedWrite(q4_x4 A, q4_x4 B, v1u_x4 Mask, v1u_x4 NotMask)
{
    q4_x4 Result = {};
    Result.x = MaskedWrite(A.x, B.x, Mask, NotMask);
    Result.y = MaskedWrite(A.y, B.y, Mask, NotMask);
    Result.z = MaskedWrite(A.z, B.z, Mask, NotMask);
    Result.z = MaskedWrite(A.w, B.w, Mask, NotMask);
    return Result;
}

inline q4_x4 MaskedWrite(q4_x4 A, q4_x4 B, v1u_x4 Mask)
{
    q4_x4 Result = MaskedWrite(A, B, Mask, ~Mask);
    return Result;
}

// =======================================================================================================================================
// NOTE: m2_x4
// =======================================================================================================================================

//
// NOTE: Init
//

inline m2_x4 M2X4(m2 A)
{
    m2_x4 Result = {};
    Result.v[0] = V2X4(A.v[0]);
    Result.v[1] = V2X4(A.v[1]);
    return Result;
}

inline m2_x4 M2X4(v2_x4 A, v2_x4 B)
{
    m2_x4 Result = {};
    Result.v[0] = A;
    Result.v[1] = B;
    return Result;
}

// TODO: Add scalar writes

//
// NOTE: Store
//

inline void StoreAligned(m2_x4 M, m2_soa Dest, u32 Id)
{
    Assert((Id % 4) == 0);
    StoreAligned(M.v[0], Dest.v[0], Id);
    StoreAligned(M.v[1], Dest.v[1], Id);
}

// TODO: Add stride and num stored
inline void StoreScalarAligned(m2_x4 M, void* Dest)
{
    // NOTE: Write out that matrices as scalar elements, instead of concatenated in a group of 4
#if MATH_X64
    __m128 r0, r1, r2, r3;
    SSE_TRANSPOSE(M.v[0].x.x, M.v[0].y.x, M.v[1].x.x, M.v[1].y.x, r0, r1, r2, r3);

    f32* CurrPtr = (f32*)Dest;    
    _mm_store_ps(CurrPtr + 0, r0);    
    _mm_store_ps(CurrPtr + 4, r1);
    _mm_store_ps(CurrPtr + 8, r2);
    _mm_store_ps(CurrPtr + 12, r3);
#elif MATH_ARM
    u8* CurrByte = (u8*)Dest;
    for (u32 ElementId = 0; ElementId < 4; ++ElementId)
    {
        f32* WriteFloat = (f32*)CurrByte;
        WriteFloat[0] = M.v[0].x.e[ElementId];
        WriteFloat[0] = M.v[0].y.e[ElementId];
        WriteFloat[0] = M.v[1].x.e[ElementId];
        WriteFloat[0] = M.v[1].y.e[ElementId];

        CurrByte += sizeof(m2);
    }
#endif
}

//
// NOTE: Add
//

inline m2_x4 operator+(m2_x4 A, m2_x4 B)
{
    m2_x4 Result = {};
    Result.v[0] = A.v[0] + B.v[0];
    Result.v[1] = A.v[1] + B.v[1];
    return Result;
}

inline m2_x4& operator+=(m2_x4& A, m2_x4 B)
{
    A = A + B;
    return A;
}

//
// NOTE: Sub
//

inline m2_x4 operator-(m2_x4 A, m2_x4 B)
{
    m2_x4 Result = {};
    Result.v[0] = A.v[0] - B.v[0];
    Result.v[1] = A.v[1] - B.v[1];
    return Result;
}

inline m2_x4& operator-=(m2_x4& A, m2_x4 B)
{
    A = A - B;
    return A;
}

//
// NOTE: Mul
//

inline m2_x4 operator*(m2_x4 A, v1_x4 B)
{
    m2_x4 Result = {};
    Result.v[0] = A.v[0] * B;
    Result.v[1] = A.v[1] * B;
    return Result;
}

inline m2_x4& operator*=(m2_x4& A, v1_x4 B)
{
    A = A * B;
    return A;
}

inline m2_x4 operator*(v1_x4 B, m2_x4 A)
{
    m2_x4 Result = {};
    Result.v[0] = A.v[0] * B;
    Result.v[1] = A.v[1] * B;
    return Result;
}

inline v2_x4 operator*(m2_x4 A, v2_x4 B)
{
    v2_x4 Result = B.x*A.v[0] + B.y*A.v[1];
    return Result;
}

inline m2_x4 operator*(m2_x4 A, m2_x4 B)
{
    m2_x4 Result = {};
    Result.v[0] = A*B.v[0];
    Result.v[1] = A*B.v[1];
    return Result;
}

inline m2_x4& operator*=(m2_x4& A, m2_x4 B)
{
    A = A * B;
    return A;
}

//
// NOTE: Mask Ops
//

inline m2_x4 MaskedWrite(m2_x4 A, m2_x4 B, v1u_x4 Mask, v1u_x4 NotMask)
{
    m2_x4 Result = {};
    Result.v[0].x = MaskedWrite(A.v[0].x, B.v[0].x, Mask, NotMask);
    Result.v[0].y = MaskedWrite(A.v[0].y, B.v[0].y, Mask, NotMask);
    Result.v[1].x = MaskedWrite(A.v[1].x, B.v[1].x, Mask, NotMask);
    Result.v[1].y = MaskedWrite(A.v[1].y, B.v[1].y, Mask, NotMask);
    return Result;
}

inline m2_x4 MaskedWrite(m2_x4 A, m2_x4 B, v1u_x4 Mask)
{
    m2_x4 Result = MaskedWrite(A, B, Mask, ~Mask);
    return Result;
}

// =======================================================================================================================================
// NOTE: m3_x4
// =======================================================================================================================================

//
// NOTE: Init
//

inline m3_x4 M3X4(m3 A)
{
    m3_x4 Result = {};
    Result.v[0] = V3X4(A.v[0]);
    Result.v[1] = V3X4(A.v[1]);
    Result.v[2] = V3X4(A.v[2]);
    return Result;
}

inline m3_x4 M3X4(v3_x4 A, v3_x4 B, v3_x4 C)
{
    m3_x4 Result = {};
    Result.v[0] = A;
    Result.v[1] = B;
    Result.v[2] = C;
    return Result;
}

// TODO: Add scalar writes

//
// NOTE: Store
//

inline void StoreAligned(m3_x4 M, m3_soa Dest, u32 Id)
{
    Assert((Id % 4) == 0);
    StoreAligned(M.v[0], Dest.v[0], Id);
    StoreAligned(M.v[1], Dest.v[1], Id);
    StoreAligned(M.v[2], Dest.v[2], Id);
}

//
// NOTE: Add
//

inline m3_x4 operator+(m3_x4 A, m3_x4 B)
{
    m3_x4 Result = {};
    Result.v[0] = A.v[0] + B.v[0];
    Result.v[1] = A.v[1] + B.v[1];
    Result.v[2] = A.v[2] + B.v[2];
    return Result;
}

inline m3_x4& operator+=(m3_x4& A, m3_x4 B)
{
    A = A + B;
    return A;
}

//
// NOTE: Sub
//

inline m3_x4 operator-(m3_x4 A, m3_x4 B)
{
    m3_x4 Result = {};
    Result.v[0] = A.v[0] - B.v[0];
    Result.v[1] = A.v[1] - B.v[1];
    Result.v[2] = A.v[2] - B.v[2];
    return Result;
}

inline m3_x4& operator-=(m3_x4& A, m3_x4 B)
{
    A = A - B;
    return A;
}

//
// NOTE: Mul
//

inline m3_x4 operator*(m3_x4 A, v1_x4 B)
{
    m3_x4 Result = {};
    Result.v[0] = A.v[0] * B;
    Result.v[1] = A.v[1] * B;
    Result.v[2] = A.v[2] * B;
    return Result;
}

inline m3_x4& operator*=(m3_x4& A, v1_x4 B)
{
    A = A * B;
    return A;
}

inline m3_x4 operator*(v1_x4 B, m3_x4 A)
{
    m3_x4 Result = {};
    Result.v[0] = A.v[0] * B;
    Result.v[1] = A.v[1] * B;
    Result.v[2] = A.v[2] * B;
    return Result;
}

inline v3_x4 operator*(m3_x4 A, v3_x4 B)
{
    v3_x4 Result = B.x*A.v[0] + B.y*A.v[1] + B.z*A.v[2];
    return Result;
}

inline m3_x4 operator*(m3_x4 A, m3_x4 B)
{
    m3_x4 Result = {};
    Result.v[0] = A*B.v[0];
    Result.v[1] = A*B.v[1];
    Result.v[2] = A*B.v[2];
    return Result;
}

inline m3_x4& operator*=(m3_x4& A, m3_x4 B)
{
    A = A * B;
    return A;
}

//
// NOTE: Mask Ops
//

inline m3_x4 MaskedWrite(m3_x4 A, m3_x4 B, v1u_x4 Mask, v1u_x4 NotMask)
{
    m3_x4 Result = {};
    Result.v[0].x = MaskedWrite(A.v[0].x, B.v[0].x, Mask, NotMask);
    Result.v[0].y = MaskedWrite(A.v[0].y, B.v[0].y, Mask, NotMask);
    Result.v[0].z = MaskedWrite(A.v[0].z, B.v[0].z, Mask, NotMask);
    Result.v[1].x = MaskedWrite(A.v[1].x, B.v[1].x, Mask, NotMask);
    Result.v[1].y = MaskedWrite(A.v[1].y, B.v[1].y, Mask, NotMask);
    Result.v[1].z = MaskedWrite(A.v[1].z, B.v[1].z, Mask, NotMask);
    Result.v[2].x = MaskedWrite(A.v[2].x, B.v[2].x, Mask, NotMask);
    Result.v[2].y = MaskedWrite(A.v[2].y, B.v[2].y, Mask, NotMask);
    Result.v[2].z = MaskedWrite(A.v[2].z, B.v[2].z, Mask, NotMask);
    return Result;
}

inline m3_x4 MaskedWrite(m3_x4 A, m3_x4 B, v1u_x4 Mask)
{
    m3_x4 Result = MaskedWrite(A, B, Mask, ~Mask);
    return Result;
}

// =======================================================================================================================================
// NOTE: m4_x4
// =======================================================================================================================================

//
// NOTE: Init
//

inline m4_x4 M4X4(m4 A)
{
    m4_x4 Result = {};
    Result.v[0] = V4X4(A.v[0]);
    Result.v[1] = V4X4(A.v[1]);
    Result.v[2] = V4X4(A.v[2]);
    Result.v[3] = V4X4(A.v[3]);
    return Result;
}

inline m4_x4 M4X4(v4_x4 A, v4_x4 B, v4_x4 C, v4_x4 D)
{
    m4_x4 Result = {};
    Result.v[0] = A;
    Result.v[1] = B;
    Result.v[2] = C;
    Result.v[3] = D;
    return Result;
}

// TODO: Add Scalar Writes

//
// NOTE: Store
//

inline void StoreAligned(m4_x4 M, m4_soa Dest, u32 Id)
{
    Assert((Id % 4) == 0);
    StoreAligned(M.v[0], Dest.v[0], Id);
    StoreAligned(M.v[1], Dest.v[1], Id);
    StoreAligned(M.v[2], Dest.v[2], Id);
    StoreAligned(M.v[3], Dest.v[3], Id);
}

inline void StoreScalarAligned(m4_x4 M, void* Dest, mm Stride, u32 NumStored = 4)
{
    /*
      NOTE: This is the initial layout

        a00_0, a00_1, a00_2, a00_3, 
        a01_0, a01_1, a01_2, a01_3, 
        a02_0, a02_1, a02_2, a02_3, 
        a03_0, a03_1, a03_2, a03_3,
        
        a10_0, a10_1, a10_2, a10_3, 
        a11_0, a11_1, a11_2, a11_3, 
        a12_0, a12_1, a12_2, a12_3, 
        a13_0, a13_1, a13_2, a13_3,
        
        a20_0, a20_1, a20_2, a20_3, 
        a21_0, a21_1, a21_2, a21_3, 
        a22_0, a22_1, a22_2, a22_3,
        a23_0, a23_1, a23_2, a23_3,
        
        a30_0, a30_1, a30_2, a30_3, 
        a31_0, a31_1, a31_2, a31_3, 
        a32_0, a32_1, a32_2, a32_3,
        a33_0, a33_1, a33_2, a33_3,

        This is the layout we want:

        a00_0, a01_0, a02_0, a03_0,
        a10_0, a11_0, a12_0, a13_0,
        a20_0, a21_0, a22_0, a23_0,
        a30_0, a31_0, a32_0, a33_0,

        a00_1, a01_1, a02_1, a03_1,
        a10_1, a11_1, a12_1, a13_1,
        a20_1, a21_1, a22_1, a23_1,
        a30_1, a31_1, a32_1, a33_1,

        a00_2, a01_2, a02_2, a03_2,
        a10_2, a11_2, a12_2, a13_2,
        a20_2, a21_2, a22_2, a23_2,
        a30_2, a31_2, a32_2, a33_2,

        a00_3, a01_3, a02_3, a03_3,
        a10_3, a11_3, a12_3, a13_3,
        a20_3, a21_3, a22_3, a23_3,
        a30_3, a31_3, a32_3, a33_3,
        
     */
    
    // NOTE: Write out that matrices as scalar elements, instead of concatenated in a group of 4
#if MATH_X64
    Assert(Stride >= sizeof(__m128)*4);
    __m128 Row0[4], Row1[4], Row2[4], Row3[4];
    SSE_TRANSPOSE(M.v[0].x.x, M.v[0].y.x, M.v[0].z.x, M.v[0].w.x, Row0[0], Row0[1], Row0[2], Row0[3]);
    SSE_TRANSPOSE(M.v[1].x.x, M.v[1].y.x, M.v[1].z.x, M.v[1].w.x, Row1[0], Row1[1], Row1[2], Row1[3]);
    SSE_TRANSPOSE(M.v[2].x.x, M.v[2].y.x, M.v[2].z.x, M.v[2].w.x, Row2[0], Row2[1], Row2[2], Row2[3]);
    SSE_TRANSPOSE(M.v[3].x.x, M.v[3].y.x, M.v[3].z.x, M.v[3].w.x, Row3[0], Row3[1], Row3[2], Row3[3]);

    u8* CurrPtr = (u8*)Dest;

    for (u32 ElementId = 0; ElementId < NumStored; ++ElementId)
    {
        _mm_store_ps((f32*)CurrPtr + 0, Row0[ElementId]);
        _mm_store_ps((f32*)CurrPtr + 4, Row1[ElementId]);
        _mm_store_ps((f32*)CurrPtr + 8, Row2[ElementId]);
        _mm_store_ps((f32*)CurrPtr + 12, Row3[ElementId]);
        CurrPtr += Stride;
    }
#elif MATH_ARM
    u8* CurrPtr = (u8*)Dest;
    for (u32 ElementId = 0; ElementId < NumStored; ++ElementId)
    {
        f32* WritePtr = (f32*)CurrPtr;
        WritePtr[0] = M.v[0].x.e[ElementId];
        WritePtr[1] = M.v[0].y.e[ElementId];
        WritePtr[2] = M.v[0].z.e[ElementId];
        WritePtr[3] = M.v[0].w.e[ElementId];

        WritePtr[4] = M.v[1].x.e[ElementId];
        WritePtr[5] = M.v[1].y.e[ElementId];
        WritePtr[6] = M.v[1].z.e[ElementId];
        WritePtr[7] = M.v[1].w.e[ElementId];

        WritePtr[8] = M.v[2].x.e[ElementId];
        WritePtr[9] = M.v[2].y.e[ElementId];
        WritePtr[10] = M.v[2].z.e[ElementId];
        WritePtr[11] = M.v[2].w.e[ElementId];

        WritePtr[12] = M.v[3].x.e[ElementId];
        WritePtr[13] = M.v[3].y.e[ElementId];
        WritePtr[14] = M.v[3].z.e[ElementId];
        WritePtr[15] = M.v[3].w.e[ElementId];

        CurrPtr += Stride;
    }
#endif
}

//
// NOTE: Add
//

inline m4_x4 operator+(m4_x4 A, m4_x4 B)
{
    m4_x4 Result = {};
    Result.v[0] = A.v[0] + B.v[0];
    Result.v[1] = A.v[1] + B.v[1];
    Result.v[2] = A.v[2] + B.v[2];
    Result.v[3] = A.v[3] + B.v[3];
    return Result;
}

inline m4_x4& operator+=(m4_x4& A, m4_x4 B)
{
    A = A + B;
    return A;
}

//
// NOTE: Sub
//

inline m4_x4 operator-(m4_x4 A, m4_x4 B)
{
    m4_x4 Result = {};
    Result.v[0] = A.v[0] - B.v[0];
    Result.v[1] = A.v[1] - B.v[1];
    Result.v[2] = A.v[2] - B.v[2];
    Result.v[3] = A.v[3] - B.v[3];
    return Result;
}

inline m4_x4& operator-=(m4_x4& A, m4_x4 B)
{
    A = A - B;
    return A;
}

//
// NOTE: Mul
//

inline m4_x4 operator*(m4_x4 A, v1_x4 B)
{
    m4_x4 Result = {};
    Result.v[0] = A.v[0] * B;
    Result.v[1] = A.v[1] * B;
    Result.v[2] = A.v[2] * B;
    Result.v[3] = A.v[3] * B;
    return Result;
}

inline m4_x4& operator*=(m4_x4& A, v1_x4 B)
{
    A = A * B;
    return A;
}

inline m4_x4 operator*(v1_x4 B, m4_x4 A)
{
    m4_x4 Result = {};
    Result.v[0] = A.v[0] * B;
    Result.v[1] = A.v[1] * B;
    Result.v[2] = A.v[2] * B;
    Result.v[3] = A.v[3] * B;
    return Result;
}

inline v4_x4 operator*(m4_x4 A, v4_x4 B)
{
    v4_x4 Result = B.x*A.v[0] + B.y*A.v[1] + B.z*A.v[2] + B.w*A.v[3];
    return Result;
}

inline m4_x4 operator*(m4_x4 A, m4_x4 B)
{
    m4_x4 Result = {};
    Result.v[0] = A*B.v[0];
    Result.v[1] = A*B.v[1];
    Result.v[2] = A*B.v[2];
    Result.v[3] = A*B.v[3];
    return Result;
}

inline m4_x4& operator*=(m4_x4& A, m4_x4 B)
{
    A = A * B;
    return A;
}

//
// NOTE: Mask Ops
//

inline m4_x4 MaskedWrite(m4_x4 A, m4_x4 B, v1u_x4 Mask, v1u_x4 NotMask)
{
    m4_x4 Result = {};
    Result.v[0].x = MaskedWrite(A.v[0].x, B.v[0].x, Mask, NotMask);
    Result.v[0].y = MaskedWrite(A.v[0].y, B.v[0].y, Mask, NotMask);
    Result.v[0].z = MaskedWrite(A.v[0].z, B.v[0].z, Mask, NotMask);
    Result.v[0].w = MaskedWrite(A.v[0].w, B.v[0].w, Mask, NotMask);
    Result.v[1].x = MaskedWrite(A.v[1].x, B.v[1].x, Mask, NotMask);
    Result.v[1].y = MaskedWrite(A.v[1].y, B.v[1].y, Mask, NotMask);
    Result.v[1].z = MaskedWrite(A.v[1].z, B.v[1].z, Mask, NotMask);
    Result.v[1].w = MaskedWrite(A.v[1].w, B.v[1].w, Mask, NotMask);
    Result.v[2].x = MaskedWrite(A.v[2].x, B.v[2].x, Mask, NotMask);
    Result.v[2].y = MaskedWrite(A.v[2].y, B.v[2].y, Mask, NotMask);
    Result.v[2].z = MaskedWrite(A.v[2].z, B.v[2].z, Mask, NotMask);
    Result.v[2].w = MaskedWrite(A.v[2].w, B.v[2].w, Mask, NotMask);
    Result.v[3].x = MaskedWrite(A.v[3].x, B.v[3].x, Mask, NotMask);
    Result.v[3].y = MaskedWrite(A.v[3].y, B.v[3].y, Mask, NotMask);
    Result.v[3].z = MaskedWrite(A.v[3].z, B.v[3].z, Mask, NotMask);
    Result.v[3].w = MaskedWrite(A.v[3].w, B.v[3].w, Mask, NotMask);
    return Result;
}

inline m4_x4 MaskedWrite(m4_x4 A, m4_x4 B, v1u_x4 Mask)
{
    m4_x4 Result = MaskedWrite(A, B, Mask, ~Mask);
    return Result;
}

// =======================================================================================================================================
// NOTE: aabb2_x4
// =======================================================================================================================================

//
// NOTE: Init
//

// TODO: Add
#if 0
inline aabb2 AabbMinMax(v2 Min, v2 Max)
{
    return { Min, Max };
}

inline aabb2 AabbMinMax(v2 Min, v2 Max)
{
    return { Min, Max };
}

inline aabb2 AabbCenterRadius(v2 Center, v2 Dim)
{
    return { Center - Dim, Center + Dim };
}
#endif

//
// NOTE: Scalar Writes
//

inline void WriteScalar(aabb2_x4* V, u32 Index, aabb2 Value)
{
    WriteScalar(&V->Min.x, Index, Value.Min.x);
    WriteScalar(&V->Min.y, Index, Value.Min.y);
    WriteScalar(&V->Max.x, Index, Value.Max.x);
    WriteScalar(&V->Max.y, Index, Value.Max.y);
}

// =======================================================================================================================================
// NOTE: aabb3_x4
// =======================================================================================================================================

// =======================================================================================================================================
// NOTE: Common Math Functions 
// =======================================================================================================================================

//
// NOTE: AndNot
//

inline v1i_x4 AndNot(v1i_x4 A, v1i_x4 B)
{
    v1i_x4 Result = {};
#if MATH_X64
    Result.x = _mm_andnot_si128(A.x, B.x);
#elif MATH_ARM
    Result = (~A) & B;
#endif
    return Result;
}

inline v1_x4 AndNot(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};
#if MATH_X64
    Result.x = _mm_andnot_ps(A.x, B.x);
#elif MATH_ARM
    Result = (~A) & B;
#endif
    return Result;
}

inline v2_x4 AndNot(v2_x4 A, v2_x4 B)
{
    v2_x4 Result = {};
    Result.x = AndNot(A.x, B.x);
    Result.y = AndNot(A.y, B.y);
    return Result;
}

inline v3_x4 AndNot(v3_x4 A, v3_x4 B)
{
    v3_x4 Result = {};
    Result.x = AndNot(A.x, B.x);
    Result.y = AndNot(A.y, B.y);
    Result.z = AndNot(A.z, B.z);
    return Result;
}

inline v4_x4 AndNot(v4_x4 A, v4_x4 B)
{
    v4_x4 Result = {};
    Result.x = AndNot(A.x, B.x);
    Result.y = AndNot(A.y, B.y);
    Result.z = AndNot(A.z, B.z);
    Result.w = AndNot(A.w, B.w);
    return Result;
}

//
// NOTE: Sign
//

inline v1_x4 Sign(v1_x4 A)
{
    v1_x4 Result = A & V1X4(0x80000000);
    return Result;
}

inline v2_x4 Sign(v2_x4 A)
{
    v2_x4 Result = {};
    Result.x = Sign(A.x);
    Result.y = Sign(A.y);
    return Result;
}

inline v3_x4 Sign(v3_x4 A)
{
    v3_x4 Result = {};
    Result.x = Sign(A.x);
    Result.y = Sign(A.y);
    Result.z = Sign(A.z);
    return Result;
}

inline v4_x4 Sign(v4_x4 A)
{
    v4_x4 Result = {};
    Result.x = Sign(A.x);
    Result.y = Sign(A.y);
    Result.z = Sign(A.z);
    Result.w = Sign(A.w);
    return Result;
}

//
// NOTE: Absolute Value
//

inline v1_x4 Abs(v1_x4 A)
{
    v1_x4 Result = A & V1X4(f32(~0x80000000));
    return Result;
}

inline v2_x4 Abs(v2_x4 A)
{
    v2_x4 Result = {};
    Result.x = Abs(A.x);
    Result.y = Abs(A.y);
    return Result;
}

inline v3_x4 Abs(v3_x4 A)
{
    v3_x4 Result = {};
    Result.x = Abs(A.x);
    Result.y = Abs(A.y);
    Result.z = Abs(A.z);
    return Result;
}

inline v4_x4 Abs(v4_x4 A)
{
    v4_x4 Result = {};
    Result.x = Abs(A.x);
    Result.y = Abs(A.y);
    Result.z = Abs(A.z);
    Result.w = Abs(A.w);
    return Result;
}

//
// NOTE: Min
//

inline v1_x4 Min(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};
#if MATH_X64
    Result.x = _mm_min_ps(A.x, B.x);
#elif MATH_ARM
    Result.e[0] = Min(A.e[0], B.e[0]);
    Result.e[1] = Min(A.e[1], B.e[1]);
    Result.e[2] = Min(A.e[2], B.e[2]);
    Result.e[3] = Min(A.e[3], B.e[3]);
#endif
    return Result;
}

inline v2_x4 Min(v2_x4 A, v2_x4 B)
{
    v2_x4 Result = {};
    Result.x = Min(A.x, B.x);
    Result.y = Min(A.y, B.y);
    return Result;
}

inline v3_x4 Min(v3_x4 A, v3_x4 B)
{
    v3_x4 Result = {};
    Result.x = Min(A.x, B.x);
    Result.y = Min(A.y, B.y);
    Result.z = Min(A.z, B.z);
    return Result;
}

inline v4_x4 Min(v4_x4 A, v4_x4 B)
{
    v4_x4 Result = {};
    Result.x = Min(A.x, B.x);
    Result.y = Min(A.y, B.y);
    Result.z = Min(A.z, B.z);
    Result.w = Min(A.w, B.w);
    return Result;
}

//
// NOTE: Max
//

inline v1_x4 Max(v1_x4 A, v1_x4 B)
{
    v1_x4 Result = {};
#if MATH_X64
    Result.x = _mm_max_ps(A.x, B.x);
#elif MATH_ARM
    Result.e[0] = Max(A.e[0], B.e[0]);
    Result.e[1] = Max(A.e[1], B.e[1]);
    Result.e[2] = Max(A.e[2], B.e[2]);
    Result.e[3] = Max(A.e[3], B.e[3]);
#endif
    return Result;
}

inline v2_x4 Max(v2_x4 A, v2_x4 B)
{
    v2_x4 Result = {};
    Result.x = Max(A.x, B.x);
    Result.y = Max(A.y, B.y);
    return Result;
}

inline v3_x4 Max(v3_x4 A, v3_x4 B)
{
    v3_x4 Result = {};
    Result.x = Max(A.x, B.x);
    Result.y = Max(A.y, B.y);
    Result.z = Max(A.z, B.z);
    return Result;
}

inline v4_x4 Max(v4_x4 A, v4_x4 B)
{
    v4_x4 Result = {};
    Result.x = Max(A.x, B.x);
    Result.y = Max(A.y, B.y);
    Result.z = Max(A.z, B.z);
    Result.w = Max(A.w, B.w);
    return Result;
}

//
// NOTE: Floor
//

// TODO: This is SSE4.1
inline v1_x4 Floor(v1_x4 A)
{
    v1_x4 Result = {};
#if MATH_X64
    Result.x = _mm_floor_ps(A.x);
#elif MATH_ARM
    Result.e[0] = FloorF32(A.e[0]);
    Result.e[1] = FloorF32(A.e[1]);
    Result.e[2] = FloorF32(A.e[2]);
    Result.e[3] = FloorF32(A.e[3]);
#endif
    return Result;
}

inline v2_x4 Floor(v2_x4 A)
{
    v2_x4 Result = {};
    Result.x = Floor(A.x);
    Result.y = Floor(A.y);
    return Result;
}

inline v3_x4 Floor(v3_x4 A)
{
    v3_x4 Result = {};
    Result.x = Floor(A.x);
    Result.y = Floor(A.y);
    Result.z = Floor(A.z);
    return Result;
}

inline v4_x4 Floor(v4_x4 A)
{
    v4_x4 Result = {};
    Result.x = Floor(A.x);
    Result.y = Floor(A.y);
    Result.z = Floor(A.z);
    Result.w = Floor(A.w);
    return Result;
}

//
// NOTE: Ceil
//

// TODO: This is SSE4.1
inline v1_x4 Ceil(v1_x4 A)
{
    v1_x4 Result = {};
#if MATH_X64
    Result.x = _mm_ceil_ps(A.x);
#elif MATH_ARM
    Result.e[0] = CeilF32(A.e[0]);
    Result.e[1] = CeilF32(A.e[1]);
    Result.e[2] = CeilF32(A.e[2]);
    Result.e[3] = CeilF32(A.e[3]);
#endif
    return Result;
}

inline v2_x4 Ceil(v2_x4 A)
{
    v2_x4 Result = {};
    Result.x = Ceil(A.x);
    Result.y = Ceil(A.y);
    return Result;
}

inline v3_x4 Ceil(v3_x4 A)
{
    v3_x4 Result = {};
    Result.x = Ceil(A.x);
    Result.y = Ceil(A.y);
    Result.z = Ceil(A.z);
    return Result;
}

inline v4_x4 Ceil(v4_x4 A)
{
    v4_x4 Result = {};
    Result.x = Ceil(A.x);
    Result.y = Ceil(A.y);
    Result.z = Ceil(A.z);
    Result.w = Ceil(A.w);
    return Result;
}

//
// NOTE: Clamp
//

inline v1_x4 Clamp(v1_x4 Val, v1_x4 MinVal, v1_x4 MaxVal)
{
    v1_x4 Result = Min(MaxVal, Max(MinVal, Val));
    return Result;
}

inline v2_x4 Clamp(v2_x4 Val, v2_x4 MinVal, v2_x4 MaxVal)
{
    v2_x4 Result = {};
    Result.x = Clamp(Val.x, MinVal.x, MaxVal.x);
    Result.y = Clamp(Val.y, MinVal.y, MaxVal.y);
    return Result;
}

inline v3_x4 Clamp(v3_x4 Val, v3_x4 MinVal, v3_x4 MaxVal)
{
    v3_x4 Result = {};
    Result.x = Clamp(Val.x, MinVal.x, MaxVal.x);
    Result.y = Clamp(Val.y, MinVal.y, MaxVal.y);
    Result.z = Clamp(Val.z, MinVal.z, MaxVal.z);
    return Result;
}

inline v4_x4 Clamp(v4_x4 Val, v4_x4 MinVal, v4_x4 MaxVal)
{
    v4_x4 Result = {};
    Result.x = Clamp(Val.x, MinVal.x, MaxVal.x);
    Result.y = Clamp(Val.y, MinVal.y, MaxVal.y);
    Result.z = Clamp(Val.z, MinVal.z, MaxVal.z);
    Result.w = Clamp(Val.w, MinVal.w, MaxVal.w);
    return Result;
}

//
// NOTE: Round
//

// TODO: This is SSE4.1
inline v1_x4 Round(v1_x4 A)
{
    v1_x4 Result = {};
#if MATH_X64
    Result.x = _mm_round_ps(A.x, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
#elif MATH_ARM
    Result.e[0] = RoundToF32(A.e[0]);
    Result.e[1] = RoundToF32(A.e[1]);
    Result.e[2] = RoundToF32(A.e[2]);
    Result.e[3] = RoundToF32(A.e[3]);
#endif
    return Result;
}

inline v2_x4 Round(v2_x4 A)
{
    v2_x4 Result = {};
    Result.x = Round(A.x);
    Result.y = Round(A.y);
    return Result;
}

inline v3_x4 Round(v3_x4 A)
{
    v3_x4 Result = {};
    Result.x = Round(A.x);
    Result.y = Round(A.y);
    Result.z = Round(A.z);
    return Result;
}

inline v4_x4 Round(v4_x4 A)
{
    v4_x4 Result = {};
    Result.x = Round(A.x);
    Result.y = Round(A.y);
    Result.z = Round(A.z);
    Result.w = Round(A.w);
    return Result;
}

//
// NOTE: Lerp
//

inline v1_x4 Lerp(v1_x4 Start, v1_x4 End, v1_x4 T)
{
    v1_x4 Result = Start*(V1X4(1.0f) - T) + End*T;
    return Result;
}

inline v1_x4 Lerp(v1_x4 Start, v1_x4 End, f32 T)
{
    v1_x4 Result = Start*(V1X4(1.0f) - T) + End*T;
    return Result;
}

inline v2_x4 Lerp(v2_x4 Start, v2_x4 End, v2_x4 T)
{
    v2_x4 Result = {};
    Result.x = Lerp(Start.x, End.x, T.x);
    Result.y = Lerp(Start.y, End.y, T.y);
    return Result;
}

inline v2_x4 Lerp(v2_x4 Start, v2_x4 End, f32 T)
{
    v2_x4 Result = {};
    Result.x = Lerp(Start.x, End.x, T);
    Result.y = Lerp(Start.y, End.y, T);
    return Result;
}

inline v3_x4 Lerp(v3_x4 Start, v3_x4 End, v3_x4 T)
{
    v3_x4 Result = {};
    Result.x = Lerp(Start.x, End.x, T.x);
    Result.y = Lerp(Start.y, End.y, T.y);
    Result.z = Lerp(Start.z, End.z, T.z);
    return Result;
}

inline v3_x4 Lerp(v3_x4 Start, v3_x4 End, f32 T)
{
    v3_x4 Result = {};
    Result.x = Lerp(Start.x, End.x, T);
    Result.y = Lerp(Start.y, End.y, T);
    Result.z = Lerp(Start.z, End.z, T);
    return Result;
}

inline v4_x4 Lerp(v4_x4 Start, v4_x4 End, v4_x4 T)
{
    v4_x4 Result = {};
    Result.x = Lerp(Start.x, End.x, T.x);
    Result.y = Lerp(Start.y, End.y, T.y);
    Result.z = Lerp(Start.z, End.z, T.z);
    Result.w = Lerp(Start.w, End.w, T.w);
    return Result;
}

inline v4_x4 Lerp(v4_x4 Start, v4_x4 End, f32 T)
{
    v4_x4 Result = {};
    Result.x = Lerp(Start.x, End.x, T);
    Result.y = Lerp(Start.y, End.y, T);
    Result.z = Lerp(Start.z, End.z, T);
    Result.w = Lerp(Start.w, End.w, T);
    return Result;
}

// TODO: Pow

//
// NOTE: Square Root
//

inline v1_x4 SquareRoot(v1_x4 A)
{
    v1_x4 Result = {};
#if MATH_X64
    Result.x = _mm_sqrt_ps(A.x);
#elif MATH_ARM
    Result.e[0] = SquareRoot(A.e[0]);
    Result.e[1] = SquareRoot(A.e[1]);
    Result.e[2] = SquareRoot(A.e[2]);
    Result.e[3] = SquareRoot(A.e[3]);
#endif
    return Result;
}

inline v2_x4 SquareRoot(v2_x4 A)
{
    v2_x4 Result = {};
    Result.x = SquareRoot(A.x);
    Result.y = SquareRoot(A.y);
    return Result;
}

inline v3_x4 SquareRoot(v3_x4 A)
{
    v3_x4 Result = {};
    Result.x = SquareRoot(A.x);
    Result.y = SquareRoot(A.y);
    Result.z = SquareRoot(A.z);
    return Result;
}

inline v4_x4 SquareRoot(v4_x4 A)
{
    v4_x4 Result = {};
    Result.x = SquareRoot(A.x);
    Result.y = SquareRoot(A.y);
    Result.z = SquareRoot(A.z);
    Result.w = SquareRoot(A.w);
    return Result;
}

//
// NOTE: Sin
//

inline v1_x4 Sin(v1_x4 x)
{
    v1_x4 Result = {};
    f32* Scalar = x.e;
    Result = V1X4(Sin(Scalar[0]), Sin(Scalar[1]), Sin(Scalar[2]), Sin(Scalar[3]));

    return Result;

    // TODO: Enable this
#if 0
    // NOTE: http://gruntthepeon.free.fr/ssemath/sse_mathfun.h
    v1_x4 Result = {};

    v1i_x4 emm0, emm2;
    v1_x4 sign_bit = x;

    x = Abs(x);
    sign_bit = Sign(sign_bit);
  
    /* scale by 4/Pi */
    Result = x * V1X4(4.0f / Pi32);

    /* store the integer part of y in mm0 */
    emm2 = V1IX4Convert(Result);
    /* j=(j+1) & (~1) (see the cephes sources) */
    emm2 = emm2 + V1IX4(1);
    emm2 = emm2 & V1IX4(~1);
    Result = V1X4Convert(emm2);

    /* get the swap sign flag */
    emm0 = emm2 & V1IX4(4);
    emm0 = emm0 << 29;
    /* get the polynom selection mask 
       there is one polynom for 0 <= x <= Pi/4
       and another one for Pi/4<x<=Pi/2

       Both branches will be computed.
    */
    emm2 = emm2 & V1IX4(2);
    emm2 = emm2 == V1IX4(0);
  
    v1_x4 swap_sign_bit = V1X4Cast(emm0);
    v1_x4 poly_mask = V1X4Cast(emm2);
    sign_bit = sign_bit ^ swap_sign_bit;
  
    /* The magic pass: "Extended precision modular arithmetic" 
       x = ((x - y * DP1) - y * DP2) - y * DP3; */
    v1_x4 xmm1 = Result*V1X4(-0.78515625f);
    v1_x4 xmm2 = Result*V1X4(-2.4187564849853515625e-4f);
    v1_x4 xmm3 = Result*V1X4(-3.77489497744594108e-8f);
    x = x + xmm1 + xmm2 + xmm3;

    /* Evaluate the first polynom  (0 <= x <= Pi/4) */
    Result = V1X4(2.443315711809948E-005f);
    v1_x4 z = x*x;

    Result = Result*z;
    Result = Result + V1X4(-1.388731625493765E-003f);
    Result = Result*z;
    Result = Result + V1X4(4.166664568298827E-002f);
    Result = Result*z;
    Result = Result*z;
    v1_x4 tmp = z*V1X4(0.5f);
    Result = Result - tmp;
    Result = Result + V1X4(1.0f);
  
    /* Evaluate the second polynom  (Pi/4 <= x <= 0) */
    v1_x4 y2 = V1X4(-1.9515295891E-4f);
    y2 = y2*z;
    y2 = y2 + V1X4(8.3321608736E-3f);
    y2 = y2*z;
    y2 = y2 + V1X4(-1.6666654611E-1f);
    y2 = y2*z;
    y2 = y2*x;
    y2 = y2 + x;

    /* select the correct result from the two polynoms */  
    xmm3 = poly_mask;
    y2 = xmm3 & y2;
    Result = AndNot(xmm3, Result);
    Result = Result + y2;
    /* update the sign */
    Result = Result ^ sign_bit;

    return Result;
#endif
}

inline v2_x4 Sin(v2_x4 Angle)
{
    v2_x4 Result = {};
    Result.x = Sin(Angle.x);
    Result.y = Sin(Angle.y);
    return Result;
}

inline v3_x4 Sin(v3_x4 Angle)
{
    v3_x4 Result = {};
    Result.x = Sin(Angle.x);
    Result.y = Sin(Angle.y);
    Result.z = Sin(Angle.z);
    return Result;
}

inline v4_x4 Sin(v4_x4 Angle)
{
    v4_x4 Result = {};
    Result.x = Sin(Angle.x);
    Result.y = Sin(Angle.y);
    Result.z = Sin(Angle.z);
    Result.w = Sin(Angle.w);
    return Result;
}

//
// NOTE: Cos
//

inline v1_x4 Cos(v1_x4 x)
{
    v1_x4 Result = {};
    f32* Scalar = x.e;
    Result = V1X4(Cos(Scalar[0]), Cos(Scalar[1]), Cos(Scalar[2]), Cos(Scalar[3]));

    return Result;

    // TODO: Add in after testing
#if 0
    // NOTE: http://gruntthepeon.free.fr/ssemath/sse_mathfun.h
    v1_x4 Result = {};
    
    v1i_x4 emm0, emm2;
  
    x = Abs(x);
  
    /* scale by 4/Pi */
    Result = x * V1X4(4.0f / Pi32);
  
    /* store the integer part of y in mm0 */
    emm2 = V1IX4Convert(Result);
    /* j=(j+1) & (~1) (see the cephes sources) */
    emm2 = emm2 + V1IX4(1);
    emm2 = emm2 & V1IX4(~1);
    Result = V1X4Convert(emm2);

    emm2 = emm2 - V1IX4(2);
  
    /* get the swap sign flag */
    // TODO: This is either a error here or in sin
    emm0 = AndNot(emm2, V1IX4(4));
    emm0 = emm0 << 29;
    /* get the polynom selection mask */
    emm2 = emm2 & V1IX4(2);
    emm2 = emm2 == V1IX4(0);
  
    v1_x4 sign_bit = V1X4Cast(emm0);
    v1_x4 poly_mask = V1X4Cast(emm2);

    /* The magic pass: "Extended precision modular arithmetic" 
       x = ((x - y * DP1) - y * DP2) - y * DP3; */
    v1_x4 xmm1 = Result * V1X4(-0.78515625f);
    v1_x4 xmm2 = Result * V1X4(-2.4187564849853515625e-4f);
    v1_x4 xmm3 = Result * V1X4(-3.77489497744594108e-8f);
    x = x + xmm1;
    x = x + xmm2;
    x = x + xmm3;
  
    /* Evaluate the first polynom  (0 <= x <= Pi/4) */
    Result = V1X4(2.443315711809948E-005f);
    v1_x4 z = x * x;

    Result = Result * z;
    Result = Result + V1X4(-1.388731625493765E-003f);
    Result = Result * z;
    Result = Result + V1X4(4.166664568298827E-002f);
    Result = Result * z;
    Result = Result * z;
    v1_x4 tmp = z * V1X4(0.5f);
    Result = Result - tmp;
    Result = Result + V1X4(1.0f);
  
    /* Evaluate the second polynom  (Pi/4 <= x <= 0) */

    v1_x4 y2 = V1X4(-1.9515295891E-4f);
    y2 = y2 * z;
    y2 = y2 + V1X4(8.3321608736E-3f);
    y2 = y2 * z;
    y2 = y2 + V1X4(-1.6666654611E-1f);
    y2 = y2 * z;
    y2 = y2 * x;
    y2 = y2 + x;

    /* select the correct result from the two polynoms */  
    xmm3 = poly_mask;
    y2 = xmm3 & y2;
    Result = AndNot(xmm3, Result);
    Result = Result + y2;
    /* update the sign */
    Result = Result ^ sign_bit;

    return Result;
#endif
}

inline v2_x4 Cos(v2_x4 Angle)
{
    v2_x4 Result = {};
    Result.x = Cos(Angle.x);
    Result.y = Cos(Angle.y);
    return Result;
}

inline v3_x4 Cos(v3_x4 Angle)
{
    v3_x4 Result = {};
    Result.x = Cos(Angle.x);
    Result.y = Cos(Angle.y);
    Result.z = Cos(Angle.z);
    return Result;
}

inline v4_x4 Cos(v4_x4 Angle)
{
    v4_x4 Result = {};
    Result.x = Cos(Angle.x);
    Result.y = Cos(Angle.y);
    Result.z = Cos(Angle.z);
    Result.w = Cos(Angle.w);
    return Result;
}

// TODO: Tan
// TODO: ArcSin
// TODO: ArcCos
// TODO: ArcTan

// =======================================================================================================================================
// NOTE: Common Vector Functions
// =======================================================================================================================================

//
// NOTE: Length Squared
//

inline v1_x4 LengthSquared(v2_x4 A)
{
    v1_x4 Result = Square(A.x) + Square(A.y);
    return Result;
}

inline v1_x4 LengthSquared(v3_x4 A)
{
    v1_x4 Result = Square(A.x) + Square(A.y) + Square(A.z);
    return Result;
}

inline v1_x4 LengthSquared(v4_x4 A)
{
    v1_x4 Result = Square(A.x) + Square(A.y) + Square(A.z) + Square(A.w);
    return Result;
}

inline v1_x4 LengthSquared(q4_x4 Q)
{
    v1_x4 Result = Square(Q.x) + Square(Q.y) + Square(Q.z) + Square(Q.w);
    return Result;
}

//
// NOTE: Length
//

inline v1_x4 Length(v2_x4 A)
{
    v1_x4 Result = SquareRoot(LengthSquared(A));
    return Result;
}

inline v1_x4 Length(v3_x4 A)
{
    v1_x4 Result = SquareRoot(LengthSquared(A));
    return Result;
}

inline v1_x4 Length(v4_x4 A)
{
    v1_x4 Result = SquareRoot(LengthSquared(A));
    return Result;
}

inline v1_x4 Length(q4_x4 Q)
{
    v1_x4 Result = SquareRoot(LengthSquared(Q));
    return Result;
}

//
// NOTE: Normalize
//

inline v2_x4 Normalize(v2_x4 A)
{
    v1_x4 VecLength = Length(A);

    f32* Scalars = VecLength.e;
    Assert(Scalars[0] != 0.0f);
    Assert(Scalars[1] != 0.0f);
    Assert(Scalars[2] != 0.0f);
    Assert(Scalars[3] != 0.0f);

    v2_x4 Result = A / VecLength;
    return Result;
}

inline v3_x4 Normalize(v3_x4 A)
{
    v1_x4 VecLength = Length(A);
    
    f32* Scalars = VecLength.e;
    Assert(Scalars[0] != 0.0f);
    Assert(Scalars[1] != 0.0f);
    Assert(Scalars[2] != 0.0f);
    Assert(Scalars[3] != 0.0f);

    v3_x4 Result = A / VecLength;
    return Result;
}

inline v4_x4 Normalize(v4_x4 A)
{
    v1_x4 VecLength = Length(A);
    
    f32* Scalars = VecLength.e;
    Assert(Scalars[0] != 0.0f);
    Assert(Scalars[1] != 0.0f);
    Assert(Scalars[2] != 0.0f);
    Assert(Scalars[3] != 0.0f);

    v4_x4 Result = A / VecLength;
    return Result;
}

inline q4_x4 Normalize(q4_x4 A)
{
    v1_x4 VecLength = Length(A);
        
    f32* Scalars = VecLength.e;
    Assert(Scalars[0] != 0.0f);
    Assert(Scalars[1] != 0.0f);
    Assert(Scalars[2] != 0.0f);
    Assert(Scalars[3] != 0.0f);

    q4_x4 Result = A / VecLength;
    return Result;
}

//
// NOTE: Normalize Safe (makes sure we return 0 if length is 0)
// TODO: Do we really want this normalize function?
//

// TODO: Add masking
#if 0
inline v2_x4 NormalizeSafe(v2_x4 A)
{
    v2_x4 Result = {};
    v1_x4 VecLength = Length(A);
    if (VecLength != 0.0f)
    {
        Result = A / VecLength;
    }
    
    return Result;
}

inline v3_x4 NormalizeSafe(v3_x4 A)
{
    v3_x4 Result = {};
    v1_x4 VecLength = Length(A);
    if (VecLength != 0.0f)
    {
        Result = A / VecLength;
    }
    
    return Result;
}

inline v4_x4 NormalizeSafe(v4_x4 A)
{
    v4_x4 Result = {};
    v1_x4 VecLength = Length(A);
    if (VecLength != 0.0f)
    {
        Result = A / VecLength;
    }
    
    return Result;
}

inline q4_x4 NormalizeSafe(q4_x4 A)
{
    q4_x4 Result = {};
    v1_x4 VecLength = Length(A);
    if (VecLength != 0.0f)
    {
        Result = A / VecLength;
    }
    
    return Result;
}
#endif

//
// NOTE: Dot
//

inline v1_x4 Dot(v2_x4 A, v2_x4 B)
{
    v1_x4 Result = A.x*B.x + A.y*B.y;
    return Result;
}

inline v1_x4 Dot(v3_x4 A, v3_x4 B)
{
    v1_x4 Result = A.x*B.x + A.y*B.y + A.z*B.z;
    return Result;
}

inline v1_x4 Dot(v4_x4 A, v4_x4 B)
{
    v1_x4 Result = A.x*B.x + A.y*B.y + A.z*B.z + A.w*B.w;
    return Result;
}

// =======================================================================================================================================
// NOTE: Common Matrix Functions
// =======================================================================================================================================

//
// NOTE: Identity
//

inline m2_x4 M2X4Identity()
{
    m2_x4 Result = M2X4(M2Identity());
    return Result;
}

inline m3_x4 M3X4Identity()
{
    m3_x4 Result = M3X4(M3Identity());
    return Result;
}

inline m4_x4 M4X4Identity()
{
    m4_x4 Result = M4X4(M4Identity());
    return Result;
}

// TODO: Transpose
// TODO: Inverse

//
// NOTE: Scale Matrix
//

inline m2_x4 M2Scale(v1_x4 X, v1_x4 Y)
{
    m2_x4 Result = M2X4(V2X4(X, V1X4(0.0f)),
                        V2X4(V1X4(0.0f), Y));
    return Result;
}

inline m2_x4 M2Scale(v2_x4 Dim)
{
    m2_x4 Result = M2Scale(Dim.x, Dim.y);
    return Result;
}

inline m3_x4 M3Scale(v1_x4 X, v1_x4 Y, v1_x4 Z)
{
    m3_x4 Result = M3X4(V3X4(X, V1X4(0.0f), V1X4(0.0f)),
                        V3X4(V1X4(0.0f), Y, V1X4(0.0f)),
                        V3X4(V1X4(0.0f), V1X4(0.0f), Z));
    return Result;
}

inline m3_x4 M3Scale(v3_x4 Dim)
{
    m3_x4 Result = M3Scale(Dim.x, Dim.y, Dim.z);
    return Result;
}

inline m4_x4 M4Scale(v1_x4 X, v1_x4 Y, v1_x4 Z, v1_x4 W)
{
    m4_x4 Result = M4X4(V4X4(X, V1X4(0.0f), V1X4(0.0f), V1X4(0.0f)),
                        V4X4(V1X4(0.0f), Y, V1X4(0.0f), V1X4(0.0f)),
                        V4X4(V1X4(0.0f), V1X4(0.0f), Z, V1X4(0.0f)),
                        V4X4(V1X4(0.0f), V1X4(0.0f), V1X4(0.0f), W));
    return Result;
}

inline m4_x4 M4Scale(v3_x4 Dim)
{
    m4_x4 Result = M4Scale(Dim.x, Dim.y, Dim.z, V1X4(1.0f));
    return Result;
}

inline m4_x4 M4Scale(v4_x4 Dim)
{
     m4_x4 Result = M4Scale(Dim.x, Dim.y, Dim.z, Dim.w);
     return Result;
}

//
// NOTE: Pos
//

inline m3_x4 M3Pos(v2_x4 Pos)
{
    m3_x4 Result = M3X4(V3X4(1, 0, 0),
                        V3X4(0, 1, 0),
                        V3X4(Pos, V1X4(1)));
    return Result;
}

inline m4_x4 M4Pos(v3_x4 Pos)
{
    m4_x4 Result = M4X4(V4X4(1, 0, 0, 0),
                        V4X4(0, 1, 0, 0),
                        V4X4(0, 0, 1, 0),
                        V4X4(Pos, V1X4(1)));
    return Result;
}

//
// NOTE: Translate
//

inline m3_x4 M3Translate(m3_x4 Mat, v2_x4 Pos)
{
    m3_x4 Result = M3X4(Mat.v[0], Mat.v[1], V3X4(Mat.v[2].xy + Pos, Mat.v[2].z));
    return Result;
}

inline m4_x4 M4Translate(m4_x4 Mat, v3_x4 Pos)
{
    m4_x4 Result = M4X4(Mat.v[0], Mat.v[1], Mat.v[2], V4X4(Mat.v[3].xyz + Pos, Mat.v[3].w));
    return Result;
}

//
// NOTE: Rotation Matrix
//

inline m3_x4 M3Rotation(v1_x4 Angle)
{
    v1_x4 CAngle = Cos(Angle);
    v1_x4 SAngle = Sin(Angle);
    m3_x4 Result = M3X4(V3X4(CAngle, SAngle, V1X4(0.0f)), V3X4(-SAngle, CAngle, V1X4(0.0f)), V3X4(0.0f, 0.0f, 1.0f));
    return Result;
}

inline m4_x4 M4Rotation(v1_x4 AngleX, v1_x4 AngleY, v1_x4 AngleZ)
{
    m4_x4 RotX = M4X4(V4X4(1.0f, 0.0f, 0.0f, 0.0f),
                      V4X4(V1X4(0.0f), Cos(AngleX), Sin(AngleX), V1X4(0.0f)),
                      V4X4(V1X4(0.0f), -Sin(AngleX), Cos(AngleX), V1X4(0.0f)),
                      V4X4(0.0f, 0.0f, 0.0f, 1.0f));
    m4_x4 RotY = M4X4(V4X4(Cos(AngleY), V1X4(0.0f), -Sin(AngleY), V1X4(0.0f)),
                      V4X4(0.0f, 1.0f, 0.0f, 0.0f),
                      V4X4(Sin(AngleY), V1X4(0.0f), Cos(AngleY), V1X4(0.0f)),
                      V4X4(0.0f, 0.0f, 0.0f, 1.0f));
    m4_x4 RotZ = M4X4(V4X4(Cos(AngleZ), Sin(AngleZ), V1X4(0.0f), V1X4(0.0f)),
                      V4X4(-Sin(AngleZ), Cos(AngleZ), V1X4(0.0f), V1X4(0.0f)),
                      V4X4(0.0f, 0.0f, 1.0f, 0.0f),
                      V4X4(0.0f, 0.0f, 0.0f, 1.0f));
    m4_x4 Result = RotX*RotY*RotZ;
    return Result;
}

inline m4_x4 RotationM4(v3_x4 Rotate)
{
    m4_x4 Result = M4Rotation(Rotate.x, Rotate.y, Rotate.z);
    return Result;
}

// =======================================================================================================================================
// NOTE: Common Aabb Functions 
// =======================================================================================================================================

//
// NOTE: Get Center
//

inline v2_x4 AabbGetCenter(aabb2_x4 A)
{
    v2_x4 Result = Lerp(A.Min, A.Max, 0.5f);
    return Result;
}

inline v3_x4 AabbGetCenter(aabb3_x4 A)
{
    v3_x4 Result = Lerp(A.Min, A.Max, 0.5f);
    return Result;
}

//
// NOTE: Get Dim
//

inline v2_x4 AabbGetDim(aabb2_x4 A)
{
    v2_x4 Result = A.Max - A.Min;
    return Result;
}

inline v3_x4 AabbGetDim(aabb3_x4 A)
{
    v3_x4 Result = A.Max - A.Min;
    return Result;
}

//
// NOTE: Get Radius
//

inline v2_x4 AabbGetRadius(aabb2_x4 A)
{
    v2_x4 Result = (A.Max - A.Min) * 0.5f;
    return Result;
}

inline v3_x4 AabbGetRadius(aabb3_x4 A)
{
    v3_x4 Result = (A.Max - A.Min) * 0.5f;
    return Result;
}

// NOTE: Soa load functions
inline v2 LoadV2(v2_soa Soa, u32 Index) { v2 Result = V2(Soa.x[Index], Soa.y[Index]); return Result; }
inline v2i LoadV2i(v2i_soa Soa, u32 Index) { v2i Result = V2i(Soa.x[Index], Soa.y[Index]); return Result; }
inline v3 LoadV3(v3_soa Soa, u32 Index) { v3 Result = V3(Soa.x[Index], Soa.y[Index], Soa.z[Index]); return Result; }
inline v4 LoadV4(v4_soa Soa, u32 Index) { v4 Result = V4(Soa.x[Index], Soa.y[Index], Soa.z[Index], Soa.w[Index]); return Result; }
inline q4 LoadQ4(q4_soa Soa, u32 Index) { q4 Result = Q4(Soa.x[Index], Soa.y[Index], Soa.z[Index], Soa.w[Index]); return Result; }
inline m2 LoadM2(m2_soa Soa, u32 Index) { m2 Result = M2(Soa.e[0][Index], Soa.e[1][Index], Soa.e[2][Index], Soa.e[3][Index]); return Result; }
inline m3 LoadM3(m3_soa Soa, u32 Index) { m3 Result = M3(Soa.e[0][Index], Soa.e[1][Index], Soa.e[2][Index], Soa.e[3][Index], Soa.e[4][Index], Soa.e[5][Index], Soa.e[6][Index], Soa.e[7][Index], Soa.e[8][Index]); return Result; }
inline m4 LoadM4(m4_soa Soa, u32 Index) { m4 Result = M4(Soa.e[0][Index], Soa.e[1][Index], Soa.e[2][Index], Soa.e[3][Index], Soa.e[4][Index], Soa.e[5][Index], Soa.e[6][Index], Soa.e[7][Index], Soa.e[8][Index], Soa.e[9][Index], Soa.e[10][Index], Soa.e[11][Index], Soa.e[12][Index], Soa.e[13][Index], Soa.e[14][Index], Soa.e[15][Index]); return Result; }
inline aabb2 LoadAabb2(aabb2_soa Soa, u32 Index) { aabb2 Result = AabbMinMax(LoadV2(Soa.Min, Index), LoadV2(Soa.Max, Index)); return Result; }
inline aabb2i LoadAabb2i(aabb2i_soa Soa, u32 Index) { aabb2i Result = AabbMinMax(LoadV2i(Soa.Min, Index), LoadV2i(Soa.Max, Index)); return Result; }
inline aabb3 LoadAabb3(aabb3_soa Soa, u32 Index) { aabb3 Result = AabbMinMax(LoadV3(Soa.Min, Index), LoadV3(Soa.Max, Index)); return Result; }

// NOTE: Soa store functions
inline void StoreV2(v2_soa Soa, u32 Index, v2 A) { Soa.x[Index] = A.x; Soa.y[Index] = A.y; }
inline void StoreV2i(v2i_soa Soa, u32 Index, v2i A) { Soa.x[Index] = A.x; Soa.y[Index] = A.y; }
inline void StoreV3(v3_soa Soa, u32 Index, v3 A) { Soa.x[Index] = A.x; Soa.y[Index] = A.y; Soa.z[Index] = A.z; }
inline void StoreV4(v4_soa Soa, u32 Index, v4 A) { Soa.x[Index] = A.x; Soa.y[Index] = A.y; Soa.z[Index] = A.z; Soa.w[Index] = A.w; }
inline void StoreQ4(q4_soa Soa, u32 Index, q4 A) { Soa.x[Index] = A.x; Soa.y[Index] = A.y; Soa.z[Index] = A.z; Soa.w[Index] = A.w; }
inline void StoreM2(m2_soa Soa, u32 Index, m2 A) { StoreV2(Soa.v[0], Index, A.v[0]); StoreV2(Soa.v[1], Index, A.v[1]); }
inline void StoreM3(m3_soa Soa, u32 Index, m3 A) { StoreV3(Soa.v[0], Index, A.v[0]); StoreV3(Soa.v[1], Index, A.v[1]); StoreV3(Soa.v[2], Index, A.v[2]); }
inline void StoreM4(m4_soa Soa, u32 Index, m4 A) { StoreV4(Soa.v[0], Index, A.v[0]); StoreV4(Soa.v[1], Index, A.v[1]); StoreV4(Soa.v[2], Index, A.v[2]); StoreV4(Soa.v[3], Index, A.v[3]); }
inline void StoreAabb2(aabb2_soa Soa, u32 Index, aabb2 A) { StoreV2(Soa.Min, Index, A.Min); StoreV2(Soa.Max, Index, A.Max); }
inline void StoreAabb2i(aabb2i_soa Soa, u32 Index, aabb2i A) { StoreV2i(Soa.Min, Index, A.Min); StoreV2i(Soa.Max, Index, A.Max); }
inline void StoreAabb3(aabb3_soa Soa, u32 Index, aabb3 A) { StoreV3(Soa.Min, Index, A.Min); StoreV3(Soa.Max, Index, A.Max); }
