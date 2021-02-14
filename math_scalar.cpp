
// =======================================================================================================================================
// NOTE: Common Math Functions 
// =======================================================================================================================================

//
// NOTE: Reinterpret Cast
//

inline u32 ReinterpretU32(f32 Val)
{
    u32 Result = reinterpret_cast<u32&>(Val);
    return Result;
}

inline u32 ReinterpretU32(i32 Val)
{
    u32 Result = reinterpret_cast<u32&>(Val);
    return Result;
}

inline i32 ReinterpretI32(f32 Val)
{
    i32 Result = reinterpret_cast<i32&>(Val);
    return Result;
}

inline i32 ReinterpretI32(u32 Val)
{
    i32 Result = reinterpret_cast<i32&>(Val);
    return Result;
}

inline f32 ReinterpretF32(u32 Val)
{
    f32 Result = reinterpret_cast<f32&>(Val);
    return Result;
}

inline f32 ReinterpretF32(i32 Val)
{
    f32 Result = reinterpret_cast<f32&>(Val);
    return Result;
}

//
// NOTE: IsNan
//

inline b32 IsNan(f32 Val)
{
    b32 Result = Val != Val;
    return Result;
}

// NOTE: Sign function
inline u8 Sign(u8 A) { return (0 < A) - (A < 0); }
inline u16 Sign(u16 A) { return (0 < A) - (A < 0); }
inline u32 Sign(u32 A) { return (0 < A) - (A < 0); }
inline u64 Sign(u64 A) { return (0 < A) - (A < 0); }
inline i8 Sign(i8 A) { return (0 < A) - (A < 0); }
inline i16 Sign(i16 A) { return (0 < A) - (A < 0); }
inline i32 Sign(i32 A) { return (0 < A) - (A < 0); }
inline i64 Sign(i64 A) { return (0 < A) - (A < 0); }
inline f32 Sign(f32 A) { return (f32)((0.0f < A) - (A < 0.0f)); }
inline f64 Sign(f64 A) { return (f64)((0.0 < A) - (A < 0.0)); }
inline v2 Sign(v2 A) { return V2(Sign(A.x), Sign(A.y)); }
inline v3 Sign(v3 A) { return V3(Sign(A.x), Sign(A.y), Sign(A.z)); }
inline v4 Sign(v4 A) { return V4(Sign(A.x), Sign(A.y), Sign(A.z), Sign(A.w)); }

// NOTE: Abs function
inline i8 Abs(i8 A) { return A < 0 ? -A : A; }
inline i16 Abs(i16 A) { return A < 0 ? -A : A; }
inline i32 Abs(i32 A) { return A < 0 ? -A : A; }
inline i64 Abs(i64 A) { return A < 0 ? -A : A; }
inline f32 Abs(f32 A) { return A < 0 ? -A : A; }
inline f64 Abs(f64 A) { return A < 0 ? -A : A; }
inline v2 Abs(v2 A) { return V2(Abs(A.x), Abs(A.y)); }
inline v3 Abs(v3 A) { return V3(Abs(A.x), Abs(A.y), Abs(A.z)); }
inline v4 Abs(v4 A) { return V4(Abs(A.x), Abs(A.y), Abs(A.z), Abs(A.w)); }

// NOTE: Min functions
inline u8 Min(u8 A, u8 B) { return A > B ? B : A; }
inline u16 Min(u16 A, u16 B) { return A > B ? B : A; }
inline u32 Min(u32 A, u32 B) { return A > B ? B : A; }
inline u64 Min(u64 A, u64 B) { return A > B ? B : A; }
inline i8 Min(i8 A, i8 B) { return A > B ? B : A; }
inline i16 Min(i16 A, i16 B) { return A > B ? B : A; }
inline i32 Min(i32 A, i32 B) { return A > B ? B : A; }
inline i64 Min(i64 A, i64 B) { return A > B ? B : A; }
inline f32 Min(f32 A, f32 B) { return A > B ? B : A; }
inline f64 Min(f64 A, f64 B) { return A > B ? B : A; }
inline v2 Min(v2 A, v2 B) { return V2(Min(A.x, B.x), Min(A.y, B.y)); }
inline v2i Min(v2i A, v2i B) { return V2i(Min(A.x, B.x), Min(A.y, B.y)); }
inline v3 Min(v3 A, v3 B) { return V3(Min(A.x, B.x), Min(A.y, B.y), Min(A.z, B.z)); }
inline v4 Min(v4 A, v4 B) { return V4(Min(A.x, B.x), Min(A.y, B.y), Min(A.z, B.z), Min(A.w, B.w)); }

// NOTE: Max functions
inline u8 Max(u8 A, u8 B) { return A > B ? A : B; }
inline u16 Max(u16 A, u16 B) { return A > B ? A : B; }
inline u32 Max(u32 A, u32 B) { return A > B ? A : B; }
inline u64 Max(u64 A, u64 B) { return A > B ? A : B; }
inline i8 Max(i8 A, i8 B) { return A > B ? A : B; }
inline i16 Max(i16 A, i16 B) { return A > B ? A : B; }
inline i32 Max(i32 A, i32 B) { return A > B ? A : B; }
inline i64 Max(i64 A, i64 B) { return A > B ? A : B; }
inline f32 Max(f32 A, f32 B) { return A > B ? A : B; }
inline f64 Max(f64 A, f64 B) { return A > B ? A : B; }
inline v2 Max(v2 A, v2 B) { return V2(Max(A.x, B.x), Max(A.y, B.y)); }
inline v2i Max(v2i A, v2i B) { return V2i(Max(A.x, B.x), Max(A.y, B.y)); }
inline v3 Max(v3 A, v3 B) { return V3(Max(A.x, B.x), Max(A.y, B.y), Max(A.z, B.z)); }
inline v4 Max(v4 A, v4 B) { return V4(Max(A.x, B.x), Max(A.y, B.y), Max(A.z, B.z), Max(A.w, B.w)); }

// NOTE: Floor functuons
inline u8 FloorU8(f32 Val) { Assert(Val >= 0.0f); return (u8)floorf(Val); }
inline u8 FloorU8(f64 Val) { Assert(Val >= 0.0f); return (u8)floor(Val); }
inline u16 FloorU16(f32 Val) { Assert(Val >= 0.0f); return (u16)floorf(Val); }
inline u16 FloorU16(f64 Val) { Assert(Val >= 0.0f); return (u16)floor(Val); }
inline u32 FloorU32(f32 Val) { Assert(Val >= 0.0f); return (u32)floorf(Val); }
inline u32 FloorU32(f64 Val) { Assert(Val >= 0.0f); return (u32)floor(Val); }
inline u64 FloorU64(f32 Val) { Assert(Val >= 0.0f); return (u64)floorf(Val); }
inline u64 FloorU64(f64 Val) { Assert(Val >= 0.0f); return (u64)floor(Val); }
inline i8 FloorI8(f32 Val) { return (i8)floorf(Val); }
inline i8 FloorI8(f64 Val) { return (i8)floor(Val); }
inline i16 FloorI16(f32 Val) { return (i16)floorf(Val); }
inline i16 FloorI16(f64 Val) { return (i16)floor(Val); }
inline i32 FloorI32(f32 Val) { return (i32)floorf(Val); }
inline i32 FloorI32(f64 Val) { return (i32)floor(Val); }
inline i64 FloorI64(f32 Val) { return (i64)floorf(Val); }
inline i64 FloorI64(f64 Val) { return (i64)floor(Val); }
inline f32 FloorF32(f32 Val) { return floorf(Val); }
inline f64 FloorF64(f64 Val) { return floor(Val); }
inline v2 FloorV2(v2 A) { return V2(FloorF32(A.x), FloorF32(A.y)); }
inline v3 FloorV3(v3 A) { return V3(FloorF32(A.x), FloorF32(A.y), FloorF32(A.z)); }
inline v4 FloorV4(v4 A) { return V4(FloorF32(A.x), FloorF32(A.y), FloorF32(A.z), FloorF32(A.w)); }

// NOTE: Ceil functuons
inline u8 CeilU8(f32 Val) { Assert(Val >= 0.0f); u8 Truncated = (u8)Val; return Truncated + (Truncated < Val); }
inline u8 CeilU8(f64 Val) { Assert(Val >= 0.0f); u8 Truncated = (u8)Val; return Truncated + (Truncated < Val); }
inline u16 CeilU16(f32 Val) { Assert(Val >= 0.0f); u16 Truncated = (u16)Val; return Truncated + (Truncated < Val); }
inline u16 CeilU16(f64 Val) { Assert(Val >= 0.0f); u16 Truncated = (u16)Val; return Truncated + (Truncated < Val); }
inline u32 CeilU32(f32 Val) { Assert(Val >= 0.0f); u32 Truncated = (u32)Val; return Truncated + (Truncated < Val); }
inline u32 CeilU32(f64 Val) { Assert(Val >= 0.0f); u32 Truncated = (u32)Val; return Truncated + (Truncated < Val); }
inline u64 CeilU64(f32 Val) { Assert(Val >= 0.0f); u64 Truncated = (u64)Val; return Truncated + (Truncated < Val); }
inline u64 CeilU64(f64 Val) { Assert(Val >= 0.0f); u64 Truncated = (u64)Val; return Truncated + (Truncated < Val); }
inline i8 CeilI8(f32 Val) { i8 Truncated = (i8)Val; return Truncated + (Truncated < Val); }
inline i8 CeilI8(f64 Val) { i8 Truncated = (i8)Val; return Truncated + (Truncated < Val); }
inline i16 CeilI16(f32 Val) { i16 Truncated = (i16)Val; return Truncated + (Truncated < Val); }
inline i16 CeilI16(f64 Val) { i16 Truncated = (i16)Val; return Truncated + (Truncated < Val); }
inline i32 CeilI32(f32 Val) { i32 Truncated = (i32)Val; return Truncated + (Truncated < Val); }
inline i32 CeilI32(f64 Val) { i32 Truncated = (i32)Val; return Truncated + (Truncated < Val); }
inline i64 CeilI64(f32 Val) { i64 Truncated = (i64)Val; return Truncated + (Truncated < Val); }
inline i64 CeilI64(f64 Val) { i64 Truncated = (i64)Val; return Truncated + (Truncated < Val); }
inline f32 CeilF32(f32 Val) { f32 Truncated = (f32)(i32)Val; return Truncated + (Truncated < Val); }
inline f64 CeilF64(f64 Val) { f64 Truncated = (f64)(i32)Val; return Truncated + (Truncated < Val); }
inline v2 CeilV2(v2 A) { return V2(CeilF32(A.x), CeilF32(A.y)); }
inline v3 CeilV3(v3 A) { return V3(CeilF32(A.x), CeilF32(A.y), CeilF32(A.z)); }
inline v4 CeilV4(v4 A) { return V4(CeilF32(A.x), CeilF32(A.y), CeilF32(A.z), CeilF32(A.w)); }

// NOTE: Clamp functions
inline u8 Clamp(u8 Val, u8 MinVal, u8 MaxVal) { return Min(MaxVal, Max(MinVal, Val)); }
inline u16 Clamp(u16 Val, u16 MinVal, u16 MaxVal) { return Min(MaxVal, Max(MinVal, Val)); }
inline u32 Clamp(u32 Val, u32 MinVal, u32 MaxVal) { return Min(MaxVal, Max(MinVal, Val)); }
inline u64 Clamp(u64 Val, u64 MinVal, u64 MaxVal) { return Min(MaxVal, Max(MinVal, Val)); }
inline i8 Clamp(i8 Val, i8 MinVal, i8 MaxVal) { return Min(MaxVal, Max(MinVal, Val)); }
inline i16 Clamp(i16 Val, i16 MinVal, i16 MaxVal) { return Min(MaxVal, Max(MinVal, Val)); }
inline i32 Clamp(i32 Val, i32 MinVal, i32 MaxVal) { return Min(MaxVal, Max(MinVal, Val)); }
inline i64 Clamp(i64 Val, i64 MinVal, i64 MaxVal) { return Min(MaxVal, Max(MinVal, Val)); }
inline f32 Clamp(f32 Val, f32 MinVal, f32 MaxVal) { return Min(MaxVal, Max(MinVal, Val)); }
inline f64 Clamp(f64 Val, f64 MinVal, f64 MaxVal) { return Min(MaxVal, Max(MinVal, Val)); }
inline v2 Clamp(v2 Val, v2 MinVal, v2 MaxVal) { return Min(MaxVal, Max(MinVal, Val)); }
inline v2i Clamp(v2i Val, v2i MinVal, v2i MaxVal) { return Min(MaxVal, Max(MinVal, Val)); }
inline v3 Clamp(v3 Val, v3 MinVal, v3 MaxVal) { return Min(MaxVal, Max(MinVal, Val)); }
inline v4 Clamp(v4 Val, v4 MinVal, v4 MaxVal) { return Min(MaxVal, Max(MinVal, Val)); }

// NOTE: Round functions
inline u8 RoundToU8(f32 A) { return (u8)(A + 0.5f); }
inline u8 RoundToU8(f64 A) { return (u8)(A + 0.5); }
inline u16 RoundToU16(f32 A) { return (u16)(A + 0.5f); }
inline u16 RoundToU16(f64 A) { return (u16)(A + 0.5); }
inline u32 RoundToU32(f32 A) { return (u32)(A + 0.5f); }
inline u32 RoundToU32(f64 A) { return (u32)(A + 0.5); }
inline u64 RoundToU64(f32 A) { return (u64)(A + 0.5f); }
inline u64 RoundToU64(f64 A) { return (u64)(A + 0.5); }
inline i8 RoundToI8(f32 A) { return (i8)(A + 0.5f); }
inline i8 RoundToI8(f64 A) { return (i8)(A + 0.5); }
inline i16 RoundToI16(f32 A) { return (i16)(A + 0.5f); }
inline i16 RoundToI16(f64 A) { return (i16)(A + 0.5); }
inline i32 RoundToI32(f32 A) { return (i32)(A + 0.5f); }
inline i32 RoundToI32(f64 A) { return (i32)(A + 0.5); }
inline i64 RoundToI64(f32 A) { return (i64)(A + 0.5f); }
inline i64 RoundToI64(f64 A) { return (i64)(A + 0.5); }
inline f32 RoundToF32(f32 A) { return (f32)((i32)(A + 0.5f)); }
inline f64 RoundToF64(f64 A) { return (f64)((i64)(A + 0.5)); }
inline v2 RoundToV2(v2 A) { return V2(RoundToF32(A.x), RoundToF32(A.y)); }
inline v3 RoundToV3(v3 A) { return V3(RoundToF32(A.x), RoundToF32(A.y), RoundToF32(A.z)); }
inline v4 RoundToV4(v4 A) { return V4(RoundToF32(A.x), RoundToF32(A.y), RoundToF32(A.z), RoundToF32(A.w)); }

// NOTE: Lerp functions
inline i32 Lerp(i32 Start, i32 End, f32 T) { return (i32)((f32)Start*(1.0f - T) + (f32)End*T); }
inline f32 Lerp(f32 Start, f32 End, f32 T) { return Start*(1.0f - T) + End*T; }
inline f64 Lerp(f64 Start, f64 End, f64 T) { return Start*(1.0 - T) + End*T; }
inline v2 Lerp(v2 Start, v2 End, v2 T) { return Start*(V2(1.0f, 1.0f) - T) + End*T; }
inline v2i Lerp(v2i Start, v2i End, v2 T) { return V2i(Lerp(V2(Start), V2(End), T)); }
inline v3 Lerp(v3 Start, v3 End, v3 T) { return Start*(V3(1.0f, 1.0f, 1.0f) - T) + End*T; }
inline v4 Lerp(v4 Start, v4 End, v4 T) { return Start*(V4(1.0f, 1.0f, 1.0f, 1.0f) - T) + End*T; }

// NOTE: Pow functions
inline f32 Pow(f32 Base, f32 Exp) { return (f32)pow(Base, Exp); }
inline f64 Pow(f64 Base, f64 Exp) { return (f64)pow(Base, Exp); }
inline v2 Pow(v2 Base, f32 Exp) { return V2(Pow(Base.x, Exp), Pow(Base.y, Exp)); }
inline v3 Pow(v3 Base, f32 Exp) { return V3(Pow(Base.x, Exp), Pow(Base.y, Exp), Pow(Base.z, Exp)); }
inline v4 Pow(v4 Base, f32 Exp) { return V4(Pow(Base.x, Exp), Pow(Base.y, Exp), Pow(Base.z, Exp), Pow(Base.w, Exp)); }

// TODO: Make these for vectors as well
inline f32 SquareRoot(f32 A)
{
    f32 Result = sqrtf(A);
    return Result;
}

//
// NOTE: Transcadental Functions
//

// NOTE: Sin
inline f32 Sin(f32 A)
{
    f32 Result = sinf(A);
    return Result;
}

inline v2 Sin(v2 A)
{
    v2 Result = V2(Sin(A.x), Sin(A.y));
    return Result;
}

inline v3 Sin(v3 A)
{
    v3 Result = V3(Sin(A.x), Sin(A.y), Sin(A.z));
    return Result;
}

inline v4 Sin(v4 A)
{
    v4 Result = V4(Sin(A.x), Sin(A.y), Sin(A.z), Sin(A.w));
    return Result;
}

// NOTE: Tan
inline f32 Tan(f32 A)
{
    f32 Result = tanf(A);
    return Result;
}

inline v2 Tan(v2 A)
{
    v2 Result = V2(Tan(A.x), Tan(A.y));
    return Result;
}

inline v3 Tan(v3 A)
{
    v3 Result = V3(Tan(A.x), Tan(A.y), Tan(A.z));
    return Result;
}

inline v4 Tan(v4 A)
{
    v4 Result = V4(Tan(A.x), Tan(A.y), Tan(A.z), Tan(A.w));
    return Result;
}

// NOTE: ArcSin
inline f32 ArcSin(f32 A)
{
    f32 Result = asinf(A);
    return Result;
}

inline v2 ArcSin(v2 A)
{
    v2 Result = V2(ArcSin(A.x), ArcSin(A.y));
    return Result;
}

inline v3 ArcSin(v3 A)
{
    v3 Result = V3(ArcSin(A.x), ArcSin(A.y), ArcSin(A.z));
    return Result;
}

inline v4 ArcSin(v4 A)
{
    v4 Result = V4(ArcSin(A.x), ArcSin(A.y), ArcSin(A.z), ArcSin(A.w));
    return Result;
}

// NOTE: Cos
inline f32 Cos(f32 A)
{
    f32 Result = cosf(A);
    return Result;
}

inline v2 Cos(v2 A)
{
    v2 Result = V2(Cos(A.x), Cos(A.y));
    return Result;
}

inline v3 Cos(v3 A)
{
    v3 Result = V3(Cos(A.x), Cos(A.y), Cos(A.z));
    return Result;
}

inline v4 Cos(v4 A)
{
    v4 Result = V4(Cos(A.x), Cos(A.y), Cos(A.z), Cos(A.w));
    return Result;
}

// NOTE: ArcCos
inline f32 ArcCos(f32 A)
{
    f32 Result = acosf(A);
    return Result;
}

inline v2 ArcCos(v2 A)
{
    v2 Result = V2(ArcCos(A.x), ArcCos(A.y));
    return Result;
}

inline v3 ArcCos(v3 A)
{
    v3 Result = V3(ArcCos(A.x), ArcCos(A.y), ArcCos(A.z));
    return Result;
}

inline v4 ArcCos(v4 A)
{
    v4 Result = V4(ArcCos(A.x), ArcCos(A.y), ArcCos(A.z), ArcCos(A.w));
    return Result;
}

// NOTE: ArcTan
inline f32 ArcTan(f32 X, f32 Y)
{
    f32 Result = atan2f(Y, X);
    return Result;
}

// NOTE: Exp
inline f32 Exp(f32 A)
{
    f32 Result = expf(A);
    return Result;
}

inline v2 Exp(v2 A)
{
    v2 Result = V2(Exp(A.x), Exp(A.y));
    return Result;
}

inline v3 Exp(v3 A)
{
    v3 Result = V3(Exp(A.x), Exp(A.y), Exp(A.z));
    return Result;
}

inline v4 Exp(v4 A)
{
    v4 Result = V4(Exp(A.x), Exp(A.y), Exp(A.z), Exp(A.w));
    return Result;
}


inline f32 MapIntoRange(f32 Val, f32 Min, f32 Max)
{
    Assert((Max - Min) != 0.0f);
    f32 Result = (Val - Min) / (Max - Min);
    
    return Result;
}

inline f32 DegreeToRadians(f32 Angle)
{
    f32 Result = Angle*Pi32/180.0f;
    return Result;
}

inline f32 RadiansToDegree(f32 Radians)
{
    f32 Result = Radians * 180.0f / Pi32;
    return Result;
}

inline b32 IsBetween(f32 Val, f32 A, f32 B)
{
    b32 Result = (A <= Val && Val <= B);
    return Result;
}

// =======================================================================================================================================
// NOTE: Common Scalar Functions
// =======================================================================================================================================

// NOTE: Divides functions
inline b32 Divides(u8 Num, u8 Denom) { u8 IntDiv = Num / Denom; f32 FltDiv = (f32)Num / (f32)Denom; return f32(IntDiv) == FltDiv; }
inline b32 Divides(u16 Num, u16 Denom) { u16 IntDiv = Num / Denom; f32 FltDiv = (f32)Num / (f32)Denom; return f32(IntDiv) == FltDiv; }
inline b32 Divides(u32 Num, u32 Denom) { u32 IntDiv = Num / Denom; f32 FltDiv = (f32)Num / (f32)Denom; return f32(IntDiv) == FltDiv; }
inline b32 Divides(u64 Num, u64 Denom) { u64 IntDiv = Num / Denom; f64 FltDiv = (f64)Num / (f32)Denom; return f64(IntDiv) == FltDiv; }
inline b32 Divides(i8 Num, i8 Denom) { i8 IntDiv = Num / Denom; f32 FltDiv = (f32)Num / (f32)Denom; return f32(IntDiv) == FltDiv; }
inline b32 Divides(i16 Num, i16 Denom) { i16 IntDiv = Num / Denom; f32 FltDiv = (f32)Num / (f32)Denom; return f32(IntDiv) == FltDiv; }
inline b32 Divides(i32 Num, i32 Denom) { i32 IntDiv = Num / Denom; f32 FltDiv = (f32)Num / (f32)Denom; return f32(IntDiv) == FltDiv; }
inline b32 Divides(i64 Num, i64 Denom) { i64 IntDiv = Num / Denom; f32 FltDiv = (f32)Num / (f32)Denom; return f32(IntDiv) == FltDiv; }

struct bit_scan_result
{
    b32 Found;
    u32 Index;
};

// TODO: There has to be a easier way to do this
inline bit_scan_result FindLeastSignificantSetBit(u32 Value)
{
    bit_scan_result Result = {};

#if COMPILER_MSVC
    Result.Found = _BitScanForward((unsigned long *)&Result.Index, Value);
#else
    for(u32 Test = 0; Test < 32; ++Test)
    {
        if(Value & (1 << Test))
        {
            Result.Index = Test;
            Result.Found = true;
            break;
        }
    }
#endif
    
    return(Result);
}

inline void GetParabolaFrom3Points(f32 x1, f32 y1, f32 x2, f32 y2, f32 x3, f32 y3, f32* A, f32* B, f32* C)
{
    // NOTE: https://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
    
    f32 Denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
    *A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / Denom;
    *B = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / Denom;
    *C = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / Denom;
}

// NOTE: Gameplay functions
inline q4 Q4FacingDirection(f32 Angle)
{
    // TODO: Remove adjustment
    f32 AdjustedAngle = Angle + (1.0f / 2.0f)*Pi32;
    q4 Result = Q4AxisAngle(V3(0, 0, 1), AdjustedAngle);

    return Result;
}

inline v4 HexToColor(u32 Hex)
{
    u8 Alpha = (Hex >> 24) & 0x8;
    u8 Red = (Hex >> 16) & 0x8;
    u8 Green = (Hex >> 8) & 0x8;
    u8 Blue = (Hex >> 0) & 0x8;

    v4 Result = (1.0f / 255.0f) * V4(f32(Red), f32(Green), f32(Blue), f32(Alpha));
    return Result;
}
