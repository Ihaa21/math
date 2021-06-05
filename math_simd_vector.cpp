
// =======================================================================================================================================
// NOTE: v2_x4
// =======================================================================================================================================

//
// NOTE: Init
//

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

inline v2_x4 V2X4(v1_x4 V)
{
    v2_x4 Result = {};
    Result.x = V;
    Result.y = V;
    return Result;
}

//
// NOTE: Memory Ops
//

inline v2_x4 V2X4LoadAligned(f32* X, f32* Y)
{
    v2_x4 Result = {};
    Result.x = V1X4LoadAligned(X);
    Result.y = V1X4LoadAligned(Y);
    return Result;
}

inline v2_x4 V2X4LoadAligned(v2_soa V, u32 Id)
{
    v2_x4 Result = V2X4LoadAligned(V.x + Id, V.y + Id);
    return Result;
}

inline v2_x4 V2X4LoadUnAligned(f32* X, f32* Y)
{
    v2_x4 Result = {};
    Result.x = V1X4LoadUnAligned(X);
    Result.y = V1X4LoadUnAligned(Y);
    return Result;
}

inline v2_x4 V2X4LoadUnAligned(v2_soa V, u32 Id)
{
    v2_x4 Result = V2X4LoadUnAligned(V.x + Id, V.y + Id);
    return Result;
}

inline v2_x4 V2X4Gather(f32* X, f32* Y, v1u_x4 Offset)
{
    v2_x4 Result = {};
    Result.x = V1X4Gather(X, Offset);
    Result.y = V1X4Gather(Y, Offset);
    return Result;
}

inline v2_x4 V2X4Gather(f32* X, f32* Y, v1i_x4 Offset)
{
    v2_x4 Result = {};
    Result.x = V1X4Gather(X, Offset);
    Result.y = V1X4Gather(Y, Offset);
    return Result;
}

inline v2_x4 V2X4Gather(f32* X, f32* Y, v1u_x4 PtrOffset, v1u_x4 Mask)
{
    v2_x4 Result = {};
    Result.x = V1X4Gather(X, PtrOffset, Mask);
    Result.y = V1X4Gather(Y, PtrOffset, Mask);
    return Result;
}

inline v2_x4 V2X4Gather(f32* X, f32* Y, v1i_x4 PtrOffset, v1u_x4 Mask)
{
    v2_x4 Result = {};
    Result.x = V1X4Gather(X, PtrOffset, Mask);
    Result.y = V1X4Gather(Y, PtrOffset, Mask);
    return Result;
}

// =======================================================================================================================================
// NOTE: Common Math Operators 
// =======================================================================================================================================

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

// =======================================================================================================================================
// NOTE: Common Math Functions 
// =======================================================================================================================================

//
// NOTE: Floor
//

inline v2_x4 Floor(v2_x4 A)
{
    v2_x4 Result = {};
    Result.x = Floor(A.x);
    Result.y = Floor(A.y);
    return Result;
}

//
// NOTE: Clamp
//

inline v2_x4 Clamp(v2_x4 Val, v2_x4 MinVal, v2_x4 MaxVal)
{
    v2_x4 Result = {};
    Result.x = Clamp(Val.x, MinVal.x, MaxVal.x);
    Result.y = Clamp(Val.y, MinVal.y, MaxVal.y);
    return Result;
}

//
// NOTE: Square Root
//

inline v2_x4 SquareRoot(v2_x4 A)
{
    v2_x4 Result = {};
    Result.x = SquareRoot(A.x);
    Result.y = SquareRoot(A.y);
    return Result;
}

//
// NOTE: Length Squared
//

inline v1_x4 LengthSquared(v2_x4 A)
{
    v1_x4 Result = Square(A.x) + Square(A.y);
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

//
// NOTE: Normalize
//

inline v2_x4 Normalize(v2_x4 A)
{
    v1_x4 VecLength = Length(A);
    v2_x4 Result = A / VecLength;
    return Result;
}
