
// =======================================================================================================================================
// NOTE: v2
// =======================================================================================================================================

//
// NOTE: Init
//

inline v2 V2(f32 Val)
{
    v2 Result = { Val, Val };
    return Result;
}

inline v2 V2(f32 X, f32 Y)
{
    v2 Result = { X, Y };
    return Result;
}

inline v2 V2(i32 X, i32 Y)
{
    v2 Result = { (f32)X, (f32)Y };
    return Result;
}

inline v2 V2(u32 X, u32 Y)
{
    v2 Result = { (f32)X, (f32)Y };
    return Result;
}

inline v2 V2(v2i V)
{
    v2 Result = { (f32)V.x, (f32)V.y };
    return Result;
}

//
// NOTE: Vector Operators
//

inline v2 operator+(f32 B, v2 A)
{
    v2 Result = V2(A.x + B, A.y + B);
    return Result;
}

inline v2 operator+(v2 A, f32 B)
{
    v2 Result = V2(A.x + B, A.y + B);
    return Result;
}

inline v2& operator+=(v2& A, f32 B)
{
    A = A + B;
    return A;
}

inline v2& operator+=(f32 A, v2& B)
{
    B = B + A;
    return B;
}

inline v2 operator-(f32 B, v2 A)
{
    v2 Result = V2(A.x - B, A.y - B);
    return Result;
}

inline v2 operator-(v2 A, f32 B)
{
    v2 Result = V2(A.x - B, A.y - B);
    return Result;
}

inline v2& operator-=(v2& A, f32 B)
{
    A = A - B; return A;
}

inline v2& operator-=(f32 A, v2& B)
{
    B = B - A;
    return B;
}

inline v2 operator*(f32 B, v2 A)
{
    v2 Result = V2(A.x * B, A.y * B);
    return Result;
}

inline v2 operator*(v2 A, f32 B)
{
    v2 Result = V2(A.x * B, A.y * B);
    return Result;
}

inline v2& operator*=(v2& A, f32 B)
{
    A = A * B;
    return A;
}

inline v2& operator*=(f32 A, v2& B)
{
    B = B * A;
    return B;
}

inline v2 operator/(v2 A, f32 B)
{
    v2 Result = V2(A.x / B, A.y / B);
    return Result;
}

inline v2 operator/(f32 A, v2 B)
{
    v2 Result = V2(A / B.x, A / B.y);
    return Result;
}

inline v2& operator/=(v2& A, f32 B)
{
    A = A / B;
    return A;
}

inline b32 operator==(v2 A, v2 B)
{
    return A.x == B.x && A.y == B.y;
}

inline v2 operator+(v2 A, v2 B)
{
    v2 Result = V2(A.x + B.x, A.y + B.y);
    return Result;
}

inline v2& operator+=(v2& A, v2 B)
{
    A = A + B;
    return A;
}

inline v2 operator-(v2 A, v2 B)
{
    v2 Result = V2(A.x - B.x, A.y - B.y);
    return Result;
}

inline v2& operator-=(v2& A, v2 B)
{
    A = A - B;
    return A;
}

inline v2 operator-(v2 A)
{
    v2 Result = V2(-A.x, -A.y);
    return Result;
}

inline v2 operator*(v2 A, v2 B)
{
    v2 Result = V2(A.x*B.x, A.y*B.y);
    return Result;
}

inline v2& operator*=(v2& A, v2 B)
{
    A = A * B;
    return A;
}

inline v2 operator/(v2 A, v2 B)
{
    v2 Result = V2(A.x/B.x, A.y/B.y);
    return Result;
}

//
// NOTE: Misc
//

inline f32 AngleBetweenVectors(v2 A, v2 B)
{
    f32 Angle = ArcCos(Dot(A, B)/(Length(A)*Length(B)));
    return Angle;
}

inline v2 GetPerp(v2 A)
{
    v2 Result = V2(A.y, -A.x);
    return Result;
}

// =======================================================================================================================================
// NOTE: v2i
// =======================================================================================================================================

inline v2i V2i(i32 Val) { return { Val, Val }; }
inline v2i V2i(i32 X, i32 Y) { return { X, Y }; }
inline v2i V2i(u32 X, u32 Y) { return { (i32)X, (i32)Y }; }
inline v2i V2i(f32 X, f32 Y) { return { (i32)X, (i32)Y }; }
inline v2i V2i(v2 V) { return { (i32)V.x, (i32)V.y }; }

inline v2i operator+(i32 B, v2i A) { return V2i(A.x + B, A.y + B); }
inline v2i operator+(v2i A, i32 B) { return V2i(A.x + B, A.y + B); }
inline v2i& operator+=(v2i& A, i32 B) { A = A + B; return A; }
inline v2i& operator+=(i32 A, v2i& B) { B = B + A; return B; }

inline v2i operator-(i32 B, v2i A) { return V2i(A.x - B, A.y - B); }
inline v2i operator-(v2i A, i32 B) { return V2i(A.x - B, A.y - B); }
inline v2i& operator-=(v2i& A, i32 B) { A = A - B; return A; }
inline v2i& operator-=(i32 A, v2i& B) { B = B - A; return B; }

inline v2i operator*(i32 B, v2i A) { return V2i(A.x * B, A.y * B); }
inline v2i operator*(v2i A, i32 B) { return V2i(A.x * B, A.y * B); }
inline v2i& operator*=(v2i& A, i32 B) { A = A * B; return A; }
inline v2i& operator*=(i32 A, v2i& B) { B = B * A; return B; }

inline v2i operator/(v2i A, i32 B) { return V2i(A.x / B, A.y / B); }
inline v2i& operator/=(v2i& A, i32 B) { A = A / B; return A; }

inline b32 operator==(v2i A, v2i B) { return A.x == B.x && A.y == B.y; }
inline v2i operator+(v2i A, v2i B) { return V2i(A.x + B.x, A.y + B.y); }
inline v2i& operator+=(v2i& A, v2i B) { A = A + B; return A; }

inline v2i operator-(v2i A, v2i B) { return V2i(A.x - B.x, A.y - B.y); }
inline v2i& operator-=(v2i& A, v2i B) { A = A - B; return A; }

inline v2i operator-(v2i A) { return V2i(-A.x, -A.y); }
inline v2i operator*(v2i A, v2i B) { return V2i(A.x*B.x, A.y*B.y); }
inline v2i& operator*=(v2i& A, v2i B) { A = A * B; return A; }
inline v2i operator/(v2i A, v2i B) { return V2i(A.x/B.x, A.y/B.y); }

// =======================================================================================================================================
// NOTE: v3
// =======================================================================================================================================

inline v3 V3(f32 Val) { return { Val, Val, Val }; }
inline v3 V3(f32 X, f32 Y, f32 Z) { return { X, Y, Z }; }
inline v3 V3(i32 X, i32 Y, i32 Z) { return { (f32)X, (f32)Y, (f32)Z }; }
inline v3 V3(u32 X, u32 Y, u32 Z) { return { (f32)X, (f32)Y, (f32)Z }; }
inline v3 V3(v2 V, f32 Z) { return V3(V.x, V.y, Z); }

inline v3 operator+(f32 B, v3 A) { return V3(A.x + B, A.y + B, A.z + B); }
inline v3 operator+(v3 A, f32 B) { return V3(A.x + B, A.y + B, A.z + B); }
inline v3& operator+=(v3& A, f32 B) { A = A + B; return A; }
inline v3& operator+=(f32 A, v3& B) { B = B + A; return B; }

inline v3 operator-(f32 B, v3 A) { return V3(A.x - B, A.y - B, A.z - B); }
inline v3 operator-(v3 A, f32 B) { return V3(A.x - B, A.y - B, A.z - B); }
inline v3& operator-=(v3& A, f32 B) { A = A - B; return A; }
inline v3& operator-=(f32 A, v3& B) { B = B - A; return B; }

inline v3 operator*(f32 B, v3 A) { return V3(A.x * B, A.y * B, A.z * B); }
inline v3 operator*(v3 A, f32 B) { return V3(A.x * B, A.y * B, A.z * B); }
inline v3& operator*=(v3& A, f32 B) { A = A * B; return A; }
inline v3& operator*=(f32 A, v3& B) { B = B * A; return B; }

inline v3 operator/(v3 A, f32 B) { return V3(A.x / B, A.y / B, A.z / B); }
inline v3& operator/=(v3& A, f32 B) { A = A / B; return A; }

inline b32 operator==(v3 A, v3 B) { return A.x == B.x && A.y == B.y && A.z == B.z; }
inline v3 operator+(v3 A, v3 B) { return V3(A.x + B.x, A.y + B.y, A.z + B.z); }
inline v3& operator+=(v3& A, v3 B) { A = A + B; return A; }

inline v3 operator-(v3 A, v3 B) { return V3(A.x - B.x, A.y - B.y, A.z - B.z); }
inline v3& operator-=(v3& A, v3 B) { A = A - B; return A; }

inline v3 operator-(v3 A) { return V3(-A.x, -A.y, -A.z); }
inline v3 operator*(v3 A, v3 B) { return V3(A.x*B.x, A.y*B.y, A.z*B.z); }
inline v3& operator*=(v3& A, v3 B) { A = A * B; return A; }
inline v3 operator/(v3 A, v3 B) { return V3(A.x/B.x, A.y/B.y, A.z/B.z); }

inline v3 Cross(v3 A, v3 B)
{
    v3 Result = {};

    Result.x = A.y*B.z - A.z*B.y;
    Result.y = A.z*B.x - A.x*B.z;
    Result.z = A.x*B.y - A.y*B.x;

    return Result;
}

inline v3 GetReflection(v3 Normal, v3 Dir)
{
    v3 Result = Dir - 2*(Dot(Dir, Normal))*Normal;
    return Result;
}

// =======================================================================================================================================
// NOTE: v4
// =======================================================================================================================================

inline v4 V4(f32 Val) { return { Val, Val, Val, Val }; }
inline v4 V4(f32 X, f32 Y, f32 Z, f32 W) { return { X, Y, Z, W }; }
inline v4 V4(i32 X, i32 Y, i32 Z, i32 W) { return { (f32)X, (f32)Y, (f32)Z, (f32)W }; }
inline v4 V4(u32 X, u32 Y, u32 Z, u32 W) { return { (f32)X, (f32)Y, (f32)Z, (f32)W }; }
inline v4 V4(v3 V, f32 W) { return V4(V.x, V.y, V.z, W); }

inline v4 operator+(f32 B, v4 A) { return V4(A.x + B, A.y + B, A.z + B, A.w + B); }
inline v4 operator+(v4 A, f32 B) { return V4(A.x + B, A.y + B, A.z + B, A.w + B); }
inline v4& operator+=(v4& A, f32 B) { A = A + B; return A; }
inline v4& operator+=(f32 A, v4& B) { B = B + A; return B; }

inline v4 operator-(f32 B, v4 A) { return V4(A.x - B, A.y - B, A.z - B, A.w - B); }
inline v4 operator-(v4 A, f32 B) { return V4(A.x - B, A.y - B, A.z - B, A.w - B); }
inline v4& operator-=(v4& A, f32 B) { A = A - B; return A; }
inline v4& operator-=(f32 A, v4& B) { B = B - A; return B; }

inline v4 operator*(f32 B, v4 A) { return V4(A.x * B, A.y * B, A.z * B, A.w * B); }
inline v4 operator*(v4 A, f32 B) { return V4(A.x * B, A.y * B, A.z * B, A.w * B); }
inline v4& operator*=(v4& A, f32 B) { A = A * B; return A; }
inline v4& operator*=(f32 A, v4& B) { B = B * A; return B; }

inline v4 operator/(v4 A, f32 B) { return V4(A.x / B, A.y / B, A.z / B, A.w / B); }
inline v4& operator/=(v4& A, f32 B) { A = A / B; return A; }

inline b32 operator==(v4 A, v4 B) { return A.x == B.x && A.y == B.y && A.z == B.z && A.w == B.w; }
inline v4 operator+(v4 A, v4 B) { return V4(A.x + B.x, A.y + B.y, A.z + B.z, A.w + B.w); }
inline v4& operator+=(v4& A, v4 B) { A = A + B; return A; }

inline v4 operator-(v4 A, v4 B) { return V4(A.x - B.x, A.y - B.y, A.z - B.z, A.w - B.w); }
inline v4& operator-=(v4& A, v4 B) { A = A - B; return A; }

inline v4 operator-(v4 A) { return V4(-A.x, -A.y, -A.z, -A.w); }
inline v4 operator*(v4 A, v4 B) { return V4(A.x*B.x, A.y*B.y, A.z*B.z, A.w*B.w); }
inline v4& operator*=(v4& A, v4 B) { A = A * B; return A; }
inline v4 operator/(v4 A, v4 B) { return V4(A.x/B.x, A.y/B.y, A.z/B.z, A.w/B.w); }

inline v4 SRGBToLinear(v4 Texel)
{
    v4 Result = (1.0f / 255.0f)*Texel;
    
    Result.r = Square(Result.r);
    Result.g = Square(Result.g);
    Result.b = Square(Result.b);

    return Result;
}

inline v4 LinearToSRGB(v4 Texel)
{
    v4 Result = Texel;
    
    Result.r = SquareRoot(Texel.r);
    Result.g = SquareRoot(Texel.g);
    Result.b = SquareRoot(Texel.b);

    Result *= 255.0f;

    return Result;
}

inline v4 PreMulAlpha(v4 Texel)
{
    v4 Result = Texel;
    Result.xyz *= Texel.a;

    return Result;
}

// =======================================================================================================================================
// NOTE: v4u
// =======================================================================================================================================

inline v4u V4u(u32 X, u32 Y, u32 Z, u32 W) { return { X, Y, Z, W }; }

// =======================================================================================================================================
// NOTE: q4
// =======================================================================================================================================

inline q4 Q4(f32 X, f32 Y, f32 Z, f32 W) { return { X, Y, Z, W }; }
inline q4 Q4(i32 X, i32 Y, i32 Z, i32 W) { return Q4((f32)X, (f32)Y, (f32)Z, (f32)W); }
inline q4 Q4(v3 V, f32 W) { return Q4(V.x, V.y, V.z, W); }

inline q4 operator*(f32 B, q4 A) { return Q4(A.x * B, A.y * B, A.z * B, A.w * B); }
inline q4 operator*(q4 A, f32 B) { return Q4(A.x * B, A.y * B, A.z * B, A.w * B); }
inline q4& operator*=(q4& A, f32 B) { A = A * B; return A; }
inline q4& operator*=(f32 A, q4& B) { B = B * A; return B; }

inline q4 operator/(q4 A, f32 B) { return Q4(A.xyz / B, A.w / B); }
inline q4& operator/=(q4& A, f32 B) { A = A / B; return A; }

inline q4 operator+(q4 A, q4 B) { return Q4(A.x + B.x, A.y + B.y, A.z + B.z, A.w + B.w); }
inline q4& operator+=(q4& A, q4 B) { A = A + B; return A; }

inline q4 operator-(q4 A, q4 B) { return Q4(A.x - B.x, A.y - B.y, A.z - B.z, A.w - B.w); }
inline q4& operator-=(q4& A, q4 B) { A = A - B; return A; }

inline q4 operator-(q4 A) { return Q4(-A.x, -A.y, -A.z, -A.w); }
// NOTE: Reference https://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/code/index.htm
inline q4 operator*(q4 A, q4 B) { return Q4(A.x * B.w + A.y * B.z - A.z * B.y + A.w * B.x,
                                            -A.x * B.z + A.y * B.w + A.z * B.x + A.w * B.y,
                                            A.x * B.y - A.y * B.x + A.z * B.w + A.w * B.z,
                                            -A.x * B.x - A.y * B.y - A.z * B.z + A.w * B.w); }
inline q4& operator*=(q4& A, q4 B) { A = A * B; return A; }

inline q4 Q4AxisAngle(v3 Axis, f32 Angle)
{
    q4 Result = {};
    Result.xyz = Axis*Sin(0.5f*Angle);
    Result.w = Cos(0.5f*Angle);

    return Result;
}

inline q4 EulerAnglesToQ4(f32 Yaw, f32 Pitch, f32 Roll)
{
    // NOTE: https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
    f32 CosYaw = Cos(Yaw*0.5f);
    f32 SinYaw = Sin(Yaw*0.5f);
    f32 CosPitch = Cos(Pitch*0.5f);
    f32 SinPitch = Sin(Pitch*0.5f);
    f32 CosRoll = Cos(Roll*0.5f);
    f32 SinRoll = Sin(Roll*0.5f);

    q4 Result = {};
    Result.x = CosRoll * SinYaw * CosPitch - SinRoll * CosYaw * SinPitch;
    Result.y = CosRoll * CosYaw * SinPitch + SinRoll * SinYaw * CosPitch;
    Result.z = SinRoll * CosYaw * CosPitch - CosRoll * SinYaw * SinPitch;
    Result.w = CosRoll * CosYaw * CosPitch + SinRoll * SinYaw * SinPitch;

    return Result;
}

inline q4 Conjugate(q4 Q)
{
    q4 Result = {};
    Result.x = -Q.x;
    Result.y = -Q.y;
    Result.z = -Q.z;
    Result.w = Q.w;

    return Result;
}

inline q4 Inverse(q4 Q)
{
    q4 Result = Conjugate(Q) / LengthSquared(Q);
    return Result;
}

inline v3 RotateVec(v3 Vec, q4 Q)
{
    v3 Result = ((Q*Q4(Vec.x, Vec.y, Vec.z, 0.0f))*Inverse(Q)).xyz;

    return Result;
}

inline v3 RotateVectorAroundAxis(v3 Vec, v3 Axis, f32 Angle)
{
    q4 Quat = Q4AxisAngle(Axis, Angle);
    v3 Result = RotateVec(Vec, Quat);
    return Result;
}

inline q4 ShortestQuaternionBetweenVectors(v3 V1, v3 V2)
{
    q4 Result = {};
    Result.xyz = Cross(V1, V2);
    Result.w = SquareRoot(LengthSquared(V1)*LengthSquared(V2)) + Dot(V1, V2);
    return Result;
}

// TODO: Implement
#if 0
inline q4 Lerp(q4 Start, q4 End, f32 t)
{
    q4 Result = ;
    return Result;
}

inline q4 Slerp(q4 Start, q4 End, f32 t)
{
    q4 Result = {};

    q4 D = 
    
    return Result;
}
#endif

inline m3 Q4ToM3(q4 Q)
{
    m3 Result = {};
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

inline m4 Q4ToM4(q4 Q)
{
    m4 Result = {};
    Result.v[0].x = 1.0f - 2.0f*Square(Q.y) - 2.0f*Square(Q.z);
    Result.v[0].y = 2.0f*Q.x*Q.y + 2.0f*Q.z*Q.w;
    Result.v[0].z = 2.0f*Q.x*Q.z - 2.0f*Q.y*Q.w;
    Result.v[0].w = 0.0f;

    Result.v[1].x = 2.0f*Q.x*Q.y - 2.0f*Q.z*Q.w;
    Result.v[1].y = 1.0f - 2.0f*Square(Q.x) - 2.0f*Square(Q.z);
    Result.v[1].z = 2.0f*Q.y*Q.z + 2.0f*Q.x*Q.w;
    Result.v[1].w = 0.0f;

    Result.v[2].x = 2.0f*Q.x*Q.z + 2.0f*Q.y*Q.w;
    Result.v[2].y = 2.0f*Q.y*Q.z - 2.0f*Q.x*Q.w;
    Result.v[2].z = 1.0f - 2.0f*Square(Q.x) - 2.0f*Square(Q.y);
    Result.v[2].w = 0.0f;

    Result.v[3].x = 0.0f;
    Result.v[3].y = 0.0f;
    Result.v[3].z = 0.0f;
    Result.v[3].w = 1.0f;

    return Result;
}

// =======================================================================================================================================
// NOTE: Common Vector Functions
// =======================================================================================================================================

// NOTE: Length Squared functions
inline f32 LengthSquared(v2 A) { return Square(A.x) + Square(A.y); }
inline f32 LengthSquared(v3 A) { return Square(A.x) + Square(A.y) + Square(A.z); }
inline f32 LengthSquared(v4 A) { return Square(A.x) + Square(A.y) + Square(A.z) + Square(A.w); }
inline f32 LengthSquared(q4 Q) { return Square(Q.x) + Square(Q.y) + Square(Q.z) + Square(Q.w); }

// NOTE: Length functions
inline f32 Length(v2 A) { return SquareRoot(LengthSquared(A)); }
inline f32 Length(v3 A) { return SquareRoot(LengthSquared(A)); }
inline f32 Length(v4 A) { return SquareRoot(LengthSquared(A)); }
inline f32 Length(q4 Q) { return SquareRoot(LengthSquared(Q)); }

// NOTE: Normalize functions
inline v2 Normalize(v2 A) { f32 VecLength = Length(A); Assert(VecLength != 0.0f); return A / VecLength; }
inline v3 Normalize(v3 A) { f32 VecLength = Length(A); Assert(VecLength != 0.0f); return A / VecLength; }
inline v4 Normalize(v4 A) { f32 VecLength = Length(A); Assert(VecLength != 0.0f); return A / VecLength; }
inline q4 Normalize(q4 A) { f32 VecLength = Length(A); Assert(VecLength != 0.0f); return A / VecLength; }

inline v2 NormalizeSafe(v2 A) { f32 VecLength = Length(A); if (VecLength == 0.0f) { return V2(0, 0); } return A / VecLength; }
inline v3 NormalizeSafe(v3 A) { f32 VecLength = Length(A); if (VecLength == 0.0f) { return V3(0, 0, 0); } return A / VecLength; }
inline v4 NormalizeSafe(v4 A) { f32 VecLength = Length(A); if (VecLength == 0.0f) { return V4(0, 0, 0, 0); } return A / VecLength; }
inline q4 NormalizeSafe(q4 A) { f32 VecLength = Length(A); if (VecLength == 0.0f) { return Q4(0, 0, 0, 0); } return A / VecLength; }

// NOTE: Dot product functions
inline f32 Dot(v2 A, v2 B) { return A.x*B.x + A.y*B.y; }
inline f32 Dot(v3 A, v3 B) { return A.x*B.x + A.y*B.y + A.z*B.z; }
inline f32 Dot(v4 A, v4 B) { return A.x*B.x + A.y*B.y + A.z*B.z + A.w*B.w; }
