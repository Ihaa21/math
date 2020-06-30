// =======================================================================================================================================
// NOTE: Scalar Helpers
// =======================================================================================================================================

// NOTE: Square
#define Square(a) ((a)*(a))

// =======================================================================================================================================
// NOTE: v2
// =======================================================================================================================================

//
// NOTE: Init
//

inline v2 V2(f32 X, f32 Y) { return { X, Y }; }
inline v2 V2(i32 X, i32 Y) { return { (f32)X, (f32)Y }; }
inline v2 V2(u32 X, u32 Y) { return { (f32)X, (f32)Y }; }
inline v2 V2(v2i V) { return { (f32)V.x, (f32)V.y }; }

// NOTE: Vector Scalar Addition
inline v2 operator+(f32 B, v2 A) { return V2(A.x + B, A.y + B); }
inline v2 operator+(v2 A, f32 B) { return V2(A.x + B, A.y + B); }
inline v2& operator+=(v2& A, f32 B) { A = A + B; return A; }
inline v2& operator+=(f32 A, v2& B) { B = B + A; return B; }

// NOTE: Vector Scalar Subtraction
inline v2 operator-(f32 B, v2 A) { return V2(A.x - B, A.y - B); }
inline v2 operator-(v2 A, f32 B) { return V2(A.x - B, A.y - B); }
inline v2& operator-=(v2& A, f32 B) { A = A - B; return A; }
inline v2& operator-=(f32 A, v2& B) { B = B - A; return B; }

// NOTE: Vector Scalar Multiplication
inline v2 operator*(f32 B, v2 A) { return V2(A.x * B, A.y * B); }
inline v2 operator*(v2 A, f32 B) { return V2(A.x * B, A.y * B); }
inline v2& operator*=(v2& A, f32 B) { A = A * B; return A; }
inline v2& operator*=(f32 A, v2& B) { B = B * A; return B; }

// NOTE: Vector Scalar Division
inline v2 operator/(v2 A, f32 B) { return V2(A.x / B, A.y / B); }
inline v2 operator/(f32 A, v2 B) { return V2(A / B.x, A / B.y); }
inline v2& operator/=(v2& A, f32 B) { A = A / B; return A; }

// NOTE: Vector Equality
inline b32 operator==(v2 A, v2 B) { return A.x == B.x && A.y == B.y; }

// NOTE: Vector Vector Addition
inline v2 operator+(v2 A, v2 B) { return V2(A.x + B.x, A.y + B.y); }
inline v2& operator+=(v2& A, v2 B) { A = A + B; return A; }

// NOTE: Vector Vector Subtraction
inline v2 operator-(v2 A, v2 B) { return V2(A.x - B.x, A.y - B.y); }
inline v2& operator-=(v2& A, v2 B) { A = A - B; return A; }

// NOTE: Vector Negation
inline v2 operator-(v2 A) { return V2(-A.x, -A.y); }
inline v2 operator*(v2 A, v2 B) { return V2(A.x*B.x, A.y*B.y); }
inline v2& operator*=(v2& A, v2 B) { A = A * B; return A; }
inline v2 operator/(v2 A, v2 B) { return V2(A.x/B.x, A.y/B.y); }

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
// NOTE: m2
// =======================================================================================================================================

// NOTE: Matrix Init
inline m2 M2(f32 C0X, f32 C0Y, f32 C1X, f32 C1Y) { return { C0X, C0Y, C1X, C1Y }; }
inline m2 M2(v2 C0, v2 C1) { m2 Result; Result.v[0] = C0; Result.v[1] = C1; return Result; }

inline m2 operator*(m2 A, f32 B) { return M2(B*A.v[0], B*A.v[1]); }
inline m2 operator*(f32 B, m2 A) { return M2(B*A.v[0], B*A.v[1]); }

// NOTE: Matrix Vec Multiplication (tested using translation mats)
inline v2 operator*(m2 B, v2 A) { return A.x*B.v[0] + A.y*B.v[1]; }
inline m2 operator*(m2 A, m2 B) { return M2(A*B.v[0], A*B.v[1]); }

// =======================================================================================================================================
// NOTE: m3
// =======================================================================================================================================

inline m3 M3(f32 C0X, f32 C0Y, f32 C0Z, f32 C1X, f32 C1Y, f32 C1Z, f32 C2X, f32 C2Y, f32 C2Z) { return { C0X, C0Y, C0Z, C1X, C1Y, C1Z, C2X, C2Y, C2Z }; }
inline m3 M3(v3 C0, v3 C1, v3 C2) { m3 Result; Result.v[0] = C0; Result.v[1] = C1; Result.v[2] = C2; return Result; }

inline m3 operator*(m3 A, f32 B) { return M3(B*A.v[0], B*A.v[1], B*A.v[2]); }
inline m3 operator*(f32 B, m3 A) { return M3(B*A.v[0], B*A.v[1], B*A.v[2]); }

inline v3 operator*(m3 B, v3 A) { return A.x*B.v[0] + A.y*B.v[1] + A.z*B.v[2]; }
inline m3 operator*(m3 A, m3 B) { return M3(A*B.v[0], A*B.v[1], A*B.v[2]); }

// =======================================================================================================================================
// NOTE: m4
// =======================================================================================================================================

inline m4 M4(f32 C0X, f32 C0Y, f32 C0Z, f32 C0W, f32 C1X, f32 C1Y, f32 C1Z, f32 C1W, f32 C2X, f32 C2Y, f32 C2Z, f32 C2W, f32 C3X, f32 C3Y, f32 C3Z, f32 C3W) { return { C0X, C0Y, C0Z, C0W, C1X, C1Y, C1Z, C1W, C2X, C2Y, C2Z, C2W, C3X, C3Y, C3Z, C3W }; }
inline m4 M4(v4 C0, v4 C1, v4 C2, v4 C3) { m4 Result; Result.v[0] = C0; Result.v[1] = C1; Result.v[2] = C2; Result.v[3] = C3; return Result; }

inline m4 operator*(m4 A, f32 B) { return M4(B*A.v[0], B*A.v[1], B*A.v[2], B*A.v[3]); }
inline m4 operator*(f32 B, m4 A) { return M4(B*A.v[0], B*A.v[1], B*A.v[2], B*A.v[3]); }

inline v4 operator*(m4 B, v4 A) { return A.x*B.v[0] + A.y*B.v[1] + A.z*B.v[2] + A.w*B.v[3]; }
inline m4 operator*(m4 A, m4 B) { return M4(A*B.v[0], A*B.v[1], A*B.v[2], A*B.v[3]); }

inline m4 LookAtM4(v3 Target, v3 Up, v3 Pos)
{
    v3 NewTarget = Normalize(Target);
    v3 NewHoriz = Normalize(Cross(Target, Up));
    v3 NewUp = Cross(NewHoriz, NewTarget);

    m4 Result = M4Identity();
    Result.v[0].xyz = NewHoriz;
    Result.v[1].xyz = NewUp;
    Result.v[2].xyz = NewTarget;
    Result.v[3].w = 1.0f;
    // TODO: Why do we transpose? Does that mean our columns are not axis vectors?
    Result = Transpose(Result);
    Result = Result*M4Pos(-Pos);
    
    return Result;
}

inline m4 PerspProjM4(f32 AspectRatio, f32 Fov, f32 Near, f32 Far)
{
    // TODO: Write out the coordinate system here + how z gets mapped
    m4 Result = {};
    Result.v[0].x = 1.0f / (AspectRatio*Tan(Fov / 2.0f));
    Result.v[1].y = -1.0f / Tan(Fov / 2.0f);
    Result.v[2].z = -Far / (Near - Far);
    Result.v[2].w = 1.0f;
    Result.v[3].z = Far*Near / (Near - Far);
    
    return Result;
}

inline m4 OrthoProjM4(f32 Left, f32 Right, f32 Top, f32 Bottom, f32 Near, f32 Far)
{
    // TODO: Write out the coordinate system here + how z gets mapped
    m4 Result = {};
    Result.v[0].x = 2.0f / (Right - Left);
    Result.v[1].y = 2.0f / (Top - Bottom);
    Result.v[2].z = -Far / (Far*(Near - Far));
    Result.v[3].x = (-Right - Left) / (Right - Left);
    Result.v[3].y = (-Top - Bottom) / (Top - Bottom);
    Result.v[3].z = Far*Near / (Far*(Near - Far));
    Result.v[3].w = 1.0f;

    return Result;
}

// =======================================================================================================================================
// NOTE: aabb2
// =======================================================================================================================================

// NOTE: Aabb functions

inline aabb2 AabbMinMax(v2 Min, v2 Max) { return { Min, Max }; }
inline aabb2 AabbCenterRadius(v2 Center, v2 Dim) { return { Center - Dim, Center + Dim }; }

// =======================================================================================================================================
// NOTE: aabb2i
// =======================================================================================================================================

inline aabb2i Aabbi(aabb2 AabbF) { return { V2i(AabbF.Min), V2i(AabbF.Max) }; }
inline aabb2i AabbiMinMax(v2i Min, v2i Max) { return { Min, Max }; }
inline aabb2i AabbiCenterRadius(v2i Center, v2i Dim) { return { Center - Dim, Center + Dim }; }

// =======================================================================================================================================
// NOTE: aabb3
// =======================================================================================================================================

inline aabb3 AabbMinMax(v3 Min, v3 Max) { return { Min, Max }; }
inline aabb3 AabbCenterRadius(v3 Center, v3 Dim) { return { Center - Dim, Center + Dim }; }

// =======================================================================================================================================
// NOTE: Plane
// =======================================================================================================================================

inline p3 P3NormalPos(v3 Normal, v3 Pos)
{
    p3 Result = {};
    Result.Normal = Normal;
    Result.d = -Dot(Normal, Pos);

    return Result;
}

inline f32 Intersect(p3 Plane, v3 Start, v3 Dir)
{
    // NOTE: https://math.stackexchange.com/questions/83990/line-and-plane-intersection-in-3d

    f32 Result = 0.0f;
    f32 ParallelCheck = Dot(Dir, Plane.Normal);
    f32 Epsilon = 0.000001f;
    if (ParallelCheck < Epsilon && ParallelCheck > -Epsilon)
    {
        // NOTE: Line is parallel to plane
        Result = NAN;
        return Result;
    }

    f32 NDotStart = Dot(Plane.Normal, Start);
    f32 NDotDir = Dot(Plane.Normal, Dir);

    return (Plane.d - NDotStart) / NDotDir;
}

// =======================================================================================================================================
// NOTE: RayCast
// =======================================================================================================================================

internal ray_cast SetupRayCast(v2 Pos, v2 DestPos)
{
    ray_cast Result = {};
    
    // NOTE: http://playtechs.blogspot.ca/2007/03/raytracing-on-grid.html
    v2 Delta = V2(Abs(DestPos.x - Pos.x), Abs(DestPos.y - Pos.y));
    Result.CurrX = i16(floor(Pos.x));
    Result.CurrY = i16(floor(Pos.y));

    Result.DeltaRecip = V2(1.0f / Delta.x, 1.0f / Delta.y);
    Result.t = 0.0f;

    Result.NumGridsToVisit = 1; // NOTE: The number of grids we will visit
    Result.IncX = 0;
    Result.IncY = 0;
    Result.NextHoriz = 0.0f;
    Result.NextVert = 0.0f;

    if (Delta.x == 0)
    {
        Result.IncX = 0;
        Result.NextHoriz = Result.DeltaRecip.x; // NOTE: Infinity
    }
    else if (DestPos.x > Pos.x)
    {
        Result.IncX = 1;
        Result.NumGridsToVisit += i16(FloorF32(DestPos.x)) - Result.CurrX;
        Result.NextHoriz = (FloorF32(Pos.x) + 1 - Pos.x) * Result.DeltaRecip.x;
    }
    else
    {
        Result.IncX = -1;
        Result.NumGridsToVisit += Result.CurrX - i16(FloorF32(DestPos.x));
        Result.NextHoriz = (Pos.x - FloorF32(Pos.x)) * Result.DeltaRecip.x;
    }

    if (Delta.y == 0)
    {
        Result.IncY = 0;
        Result.NextVert = Result.DeltaRecip.y; // NOTE: Infinity
    }
    else if (DestPos.y > Pos.y)
    {
        Result.IncY = 1;
        Result.NumGridsToVisit += i16(FloorF32(DestPos.y)) - Result.CurrY;
        Result.NextVert = (FloorF32(Pos.y) + 1 - Pos.y) * Result.DeltaRecip.y;
    }
    else
    {
        Result.IncY = -1;
        Result.NumGridsToVisit += Result.CurrY - i16(FloorF32(DestPos.y));
        Result.NextVert = (Pos.y - FloorF32(Pos.y)) * Result.DeltaRecip.y;
    }

    return Result;
}

inline ray_cast IncrementRayCast(ray_cast Ray)
{
    ray_cast Result = Ray;
    if (Result.NextVert < Result.NextHoriz)
    {
        Result.CurrY += Result.IncY;
        Result.t = Result.NextVert;
        Result.NextVert += Result.DeltaRecip.y;
    }
    else
    {
        Result.CurrX += Result.IncX;
        Result.t = Result.NextHoriz;
        Result.NextHoriz += Result.DeltaRecip.x;
    }

    return Result;
}

inline f32 Intersect(v2 RayPos, v2 RayDir, aabb2 Aabb)
{
    v2 InvRayDir = 1.0f / RayDir;
    
    f32 Tx1 = (Aabb.Min.x - RayPos.x)*InvRayDir.x;
    f32 Tx2 = (Aabb.Max.x - RayPos.x)*InvRayDir.x;
 
    f32 MinT = Min(Tx1, Tx2);
    f32 MaxT = Max(Tx1, Tx2);
 
    f32 Ty1 = (Aabb.Min.y - RayPos.y)*InvRayDir.y;
    f32 Ty2 = (Aabb.Max.y - RayPos.y)*InvRayDir.y;
 
    MinT = Max(MinT, Min(Ty1, Ty2));
    MaxT = Min(MaxT, Max(Ty1, Ty2));

    f32 ResultT = -1.0f;
    if (MaxT >= MinT)
    {
        ResultT = MinT;
    }
    
    return ResultT;
}

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
inline v2 Lerp(v2 Start, v2 End, f32 T) { return Start*(V2(1.0f, 1.0f) - T) + End*T; }
inline v2i Lerp(v2i Start, v2i End, f32 T) { return V2i(Lerp(V2(Start), V2(End), T)); }
inline v3 Lerp(v3 Start, v3 End, f32 T) { return Start*(V3(1.0f, 1.0f, 1.0f) - T) + End*T; }
inline v4 Lerp(v4 Start, v4 End, f32 T) { return Start*(V4(1.0f, 1.0f, 1.0f, 1.0f) - T) + End*T; }

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

inline f32 Sin(f32 A)
{
    f32 Result = sinf(A);
    return Result;
}

inline f32 Tan(f32 A)
{
    f32 Result = tanf(A);
    return Result;
}

inline f32 ArcSin(f32 A)
{
    f32 Result = asinf(A);
    return Result;
}

inline f32 Cos(f32 A)
{
    f32 Result = cosf(A);
    return Result;
}

inline f32 ArcCos(f32 A)
{
    f32 Result = acosf(A);
    return Result;
}

inline f32 ArcTan(f32 X, f32 Y)
{
    f32 Result = atan2f(Y, X);
    return Result;
}

inline f32 MapIntoRange(f32 Val, f32 Min, f32 Max)
{
    Assert((Max - Min) != 0.0f);
    f32 Result = (Val - Min) / (Max - Min);
    
    return Result;
}

inline f32 DegreeToRad(f32 Angle)
{
    f32 Result = Angle*Pi32/180.0f;
    return Result;
}

inline f32 RadToDegree(f32 Radians)
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
inline b32 Divides(u8 Num, u8 Denom) { u8 IntDiv = Num / Denom; f32 FltDiv = (f32)Num / (f32)Denom; return IntDiv == RoundToU8(FltDiv); }
inline b32 Divides(u16 Num, u16 Denom) { u16 IntDiv = Num / Denom; f32 FltDiv = (f32)Num / (f32)Denom; return IntDiv == RoundToU16(FltDiv); }
inline b32 Divides(u32 Num, u32 Denom) { u32 IntDiv = Num / Denom; f32 FltDiv = (f32)Num / (f32)Denom; return IntDiv == RoundToU32(FltDiv); }
inline b32 Divides(u64 Num, u64 Denom) { u64 IntDiv = Num / Denom; f32 FltDiv = (f32)Num / (f32)Denom; return IntDiv == RoundToU64(FltDiv); }
inline b32 Divides(i8 Num, i8 Denom) { i8 IntDiv = Num / Denom; f32 FltDiv = (f32)Num / (f32)Denom; return IntDiv == RoundToI8(FltDiv); }
inline b32 Divides(i16 Num, i16 Denom) { i16 IntDiv = Num / Denom; f32 FltDiv = (f32)Num / (f32)Denom; return IntDiv == RoundToI16(FltDiv); }
inline b32 Divides(i32 Num, i32 Denom) { i32 IntDiv = Num / Denom; f32 FltDiv = (f32)Num / (f32)Denom; return IntDiv == RoundToI32(FltDiv); }
inline b32 Divides(i64 Num, i64 Denom) { i64 IntDiv = Num / Denom; f32 FltDiv = (f32)Num / (f32)Denom; return IntDiv == RoundToI64(FltDiv); }

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

// =======================================================================================================================================
// NOTE: Common Matrix Functions
// =======================================================================================================================================

// NOTE: Matrix Identity
inline m2 M2Identity() { return M2(V2(1, 0), V2(0, 1)); }
inline m3 M3Identity() { return M3(V3(1, 0, 0), V3(0, 1, 0), V3(0, 0, 1)); }
inline m4 M4Identity() { return M4(V4(1, 0, 0, 0), V4(0, 1, 0, 0), V4(0, 0, 1, 0), V4(0, 0, 0, 1)); }

// NOTE: Transpose Matrix
inline m2 Transpose(m2 M) { return M2(V2(M.v[0].x, M.v[1].x),
                                      V2(M.v[0].y, M.v[1].y)); }
inline m3 Transpose(m3 M) { return M3(V3(M.v[0].x, M.v[1].x, M.v[2].x),
                                      V3(M.v[0].y, M.v[1].y, M.v[2].y),
                                      V3(M.v[0].z, M.v[1].z, M.v[2].z)); }
inline m4 Transpose(m4 M) { return M4(V4(M.v[0].x, M.v[1].x, M.v[2].x, M.v[3].x),
                                      V4(M.v[0].y, M.v[1].y, M.v[2].y, M.v[3].y),
                                      V4(M.v[0].z, M.v[1].z, M.v[2].z, M.v[3].z),
                                      V4(M.v[0].w, M.v[1].w, M.v[2].w, M.v[3].w)); }

// NOTE: Matrix Minor
inline m2 MatrixMinor(m3 M, u32 ChosenX, u32 ChosenY)
{
    m2 Result = {};

    u32 CurrElement = 0;
    for (u32 Y = 0; Y < 3; ++Y)
    {
        if (Y != ChosenY)
        {
            for (u32 X = 0; X < 3; ++X)
            {
                if (X != ChosenX)
                {
                    Result.e[CurrElement++] = M.e[Y*3 + X];
                }
            }
        }
    }
    
    return Result;
}

inline m3 MatrixMinor(m4 M, u32 ChosenX, u32 ChosenY)
{
    m3 Result = {};

    u32 CurrElement = 0;
    for (u32 Y = 0; Y < 4; ++Y)
    {
        if (Y != ChosenY)
        {
            for (u32 X = 0; X < 4; ++X)
            {
                if (X != ChosenX)
                {
                    Result.e[CurrElement++] = M.e[Y*4 + X];
                }
            }
        }
    }
    
    return Result;
}

// NOTE: Matrix cofactors
inline m2 MatrixCoFactor(m2 M)
{
    m2 Result = M;
    Result.e[0] *= -1.0f;
    Result.e[2] *= -1.0f;

    return Result;
}

inline m3 MatrixCoFactor(m3 M)
{
    m3 Result = M;
    Result.e[1] *= -1.0f;
    Result.e[3] *= -1.0f;
    Result.e[5] *= -1.0f;
    Result.e[7] *= -1.0f;

    return Result;
}

inline m4 MatrixCoFactor(m4 M)
{
    // TODO: Implement
    m4 Result = M;
    Result.e[1] *= -1.0f;
    Result.e[3] *= -1.0f;
    Result.e[5] *= -1.0f;
    Result.e[7] *= -1.0f;

    return Result;
}

// NOTE: Determinant of Matrix
inline f32 Determinant(m2 M) { return (M.e[0]*M.e[2]) - (M.e[1]*M.e[3]); }
inline f32 Determinant(m3 M) { return (M.e[0]*Determinant(M2(M.v[1].yz, M.v[2].yz))
                                       - M.e[3]*Determinant(M2(M.v[0].yz, M.v[2].yz))
                                       + M.e[6]*Determinant(M2(M.v[0].yz, M.v[1].yz))); }
inline f32 Determinant(m4 M) { return (M.e[0]*Determinant(M3(M.v[1].yzw, M.v[2].yzw, M.v[3].yzw))
                                       - M.e[4]*Determinant(M3(M.v[0].yzw, M.v[2].yzw, M.v[3].yzw))
                                       + M.e[8]*Determinant(M3(M.v[0].yzw, M.v[1].yzw, M.v[3].yzw))
                                       - M.e[12]*Determinant(M3(M.v[0].yzw, M.v[1].yzw, M.v[2].yzw))); }

// NOTE: Inverse matrix
inline m2 Inverse(m2 A)
{
    f32 Det = Determinant(A);
    Assert(Det != 0.0f);
    return (1.0f / Det)*M2(A.e[3], -A.e[2], -A.e[1], A.e[0]);
}

inline m3 Inverse(m3 A)
{
    f32 Det = Determinant(A);
    Assert(Det != 0.0f);
    return (1.0f / Det)*Transpose(M3(Determinant(M2(A.v[1].yz, A.v[2].yz)),
                                     -Determinant(M2(A.v[0].yz, A.v[2].yz)),
                                     Determinant(M2(A.v[0].yz, A.v[1].yz)),
                                     -Determinant(M2(A.e[3], A.e[5], A.e[6], A.e[8])),
                                     Determinant(M2(A.e[0], A.e[1], A.e[6], A.e[8])),
                                     -Determinant(M2(A.e[0], A.e[1], A.e[3], A.e[5])),
                                     Determinant(M2(A.v[1].xy, A.v[2].xy)),
                                     -Determinant(M2(A.v[0].xy, A.v[2].xy)),
                                     Determinant(M2(A.v[0].xy, A.v[1].xy))));
}

#if 0
inline m4 Inverse(m4 A)
{
    f32 Det = Determinant(A);
    Assert(Det != 0.0f);
    return (1.0f / Det)*Transpose(M4(Determinant(M3(A.v[1].yzw, A.v[2].yzw, A.v[3].yzw)),
                                     -Determinant(M3(A.v[0].yzw, A.v[2].yzw, A.v[3].yzw)),
                                     Determinant(M3(A.v[0].yzw, A.v[1].yzw, A.v[2].yzw)),

                                     -Determinant(M3(A.v[1].x, )),
                                     Determinant(M2(A.e[0], A.e[1], A.e[6], A.e[8])),
                                     -Determinant(M2(A.e[0], A.e[1], A.e[3], A.e[5])),

                                     Determinant(M2(A.v[1].xy, A.v[2].xy)),
                                     -Determinant(M2(A.v[0].xy, A.v[2].xy)),
                                     Determinant(M2(A.v[0].xy, A.v[1].xy))));
}
#endif

inline m4 Inverse(m4 A)
{
    m4 Result = {};
    
    f32 Inv[16];

    Inv[0] = A.e[5]  * A.e[10] * A.e[15] - 
             A.e[5]  * A.e[11] * A.e[14] - 
             A.e[9]  * A.e[6]  * A.e[15] + 
             A.e[9]  * A.e[7]  * A.e[14] +
             A.e[13] * A.e[6]  * A.e[11] - 
             A.e[13] * A.e[7]  * A.e[10];

    Inv[4] = -A.e[4]  * A.e[10] * A.e[15] + 
              A.e[4]  * A.e[11] * A.e[14] + 
              A.e[8]  * A.e[6]  * A.e[15] - 
              A.e[8]  * A.e[7]  * A.e[14] - 
              A.e[12] * A.e[6]  * A.e[11] + 
              A.e[12] * A.e[7]  * A.e[10];

    Inv[8] = A.e[4]  * A.e[9] * A.e[15] - 
             A.e[4]  * A.e[11] * A.e[13] - 
             A.e[8]  * A.e[5] * A.e[15] + 
             A.e[8]  * A.e[7] * A.e[13] + 
             A.e[12] * A.e[5] * A.e[11] - 
             A.e[12] * A.e[7] * A.e[9];

    Inv[12] = -A.e[4]  * A.e[9] * A.e[14] + 
               A.e[4]  * A.e[10] * A.e[13] +
               A.e[8]  * A.e[5] * A.e[14] - 
               A.e[8]  * A.e[6] * A.e[13] - 
               A.e[12] * A.e[5] * A.e[10] + 
               A.e[12] * A.e[6] * A.e[9];

    Inv[1] = -A.e[1]  * A.e[10] * A.e[15] + 
              A.e[1]  * A.e[11] * A.e[14] + 
              A.e[9]  * A.e[2] * A.e[15] - 
              A.e[9]  * A.e[3] * A.e[14] - 
              A.e[13] * A.e[2] * A.e[11] + 
              A.e[13] * A.e[3] * A.e[10];

    Inv[5] = A.e[0]  * A.e[10] * A.e[15] - 
             A.e[0]  * A.e[11] * A.e[14] - 
             A.e[8]  * A.e[2] * A.e[15] + 
             A.e[8]  * A.e[3] * A.e[14] + 
             A.e[12] * A.e[2] * A.e[11] - 
             A.e[12] * A.e[3] * A.e[10];

    Inv[9] = -A.e[0]  * A.e[9] * A.e[15] + 
              A.e[0]  * A.e[11] * A.e[13] + 
              A.e[8]  * A.e[1] * A.e[15] - 
              A.e[8]  * A.e[3] * A.e[13] - 
              A.e[12] * A.e[1] * A.e[11] + 
              A.e[12] * A.e[3] * A.e[9];

    Inv[13] = A.e[0]  * A.e[9] * A.e[14] - 
              A.e[0]  * A.e[10] * A.e[13] - 
              A.e[8]  * A.e[1] * A.e[14] + 
              A.e[8]  * A.e[2] * A.e[13] + 
              A.e[12] * A.e[1] * A.e[10] - 
              A.e[12] * A.e[2] * A.e[9];

    Inv[2] = A.e[1]  * A.e[6] * A.e[15] - 
             A.e[1]  * A.e[7] * A.e[14] - 
             A.e[5]  * A.e[2] * A.e[15] + 
             A.e[5]  * A.e[3] * A.e[14] + 
             A.e[13] * A.e[2] * A.e[7] - 
             A.e[13] * A.e[3] * A.e[6];

    Inv[6] = -A.e[0]  * A.e[6] * A.e[15] + 
              A.e[0]  * A.e[7] * A.e[14] + 
              A.e[4]  * A.e[2] * A.e[15] - 
              A.e[4]  * A.e[3] * A.e[14] - 
              A.e[12] * A.e[2] * A.e[7] + 
              A.e[12] * A.e[3] * A.e[6];

    Inv[10] = A.e[0]  * A.e[5] * A.e[15] - 
              A.e[0]  * A.e[7] * A.e[13] - 
              A.e[4]  * A.e[1] * A.e[15] + 
              A.e[4]  * A.e[3] * A.e[13] + 
              A.e[12] * A.e[1] * A.e[7] - 
              A.e[12] * A.e[3] * A.e[5];

    Inv[14] = -A.e[0]  * A.e[5] * A.e[14] + 
               A.e[0]  * A.e[6] * A.e[13] + 
               A.e[4]  * A.e[1] * A.e[14] - 
               A.e[4]  * A.e[2] * A.e[13] - 
               A.e[12] * A.e[1] * A.e[6] + 
               A.e[12] * A.e[2] * A.e[5];

    Inv[3] = -A.e[1] * A.e[6] * A.e[11] + 
              A.e[1] * A.e[7] * A.e[10] + 
              A.e[5] * A.e[2] * A.e[11] - 
              A.e[5] * A.e[3] * A.e[10] - 
              A.e[9] * A.e[2] * A.e[7] + 
              A.e[9] * A.e[3] * A.e[6];

    Inv[7] = A.e[0] * A.e[6] * A.e[11] - 
             A.e[0] * A.e[7] * A.e[10] - 
             A.e[4] * A.e[2] * A.e[11] + 
             A.e[4] * A.e[3] * A.e[10] + 
             A.e[8] * A.e[2] * A.e[7] - 
             A.e[8] * A.e[3] * A.e[6];

    Inv[11] = -A.e[0] * A.e[5] * A.e[11] + 
               A.e[0] * A.e[7] * A.e[9] + 
               A.e[4] * A.e[1] * A.e[11] - 
               A.e[4] * A.e[3] * A.e[9] - 
               A.e[8] * A.e[1] * A.e[7] + 
               A.e[8] * A.e[3] * A.e[5];

    Inv[15] = A.e[0] * A.e[5] * A.e[10] - 
              A.e[0] * A.e[6] * A.e[9] - 
              A.e[4] * A.e[1] * A.e[10] + 
              A.e[4] * A.e[2] * A.e[9] + 
              A.e[8] * A.e[1] * A.e[6] - 
              A.e[8] * A.e[2] * A.e[5];

    f32 Determinant = A.e[0]*Inv[0] + A.e[1]*Inv[4] + A.e[2]*Inv[8] + A.e[3]*Inv[12];
    Assert(Determinant != 0.0f);
    Determinant = 1.0f / Determinant;

    for (u32 ElementIndex = 0; ElementIndex < 16; ++ElementIndex)
    {
        Result.e[ElementIndex] = Inv[ElementIndex] * Determinant;
    }
    
    return Result;
}

// NOTE: Scale Matrix
inline m2 M2Scale(f32 X, f32 Y) { return M2(V2(X, 0.0f),
                                            V2(0.0f, Y)); }
inline m2 M2Scale(v2 Dim) { return M2Scale(Dim.x, Dim.y); }
inline m3 M3Scale(f32 X, f32 Y, f32 Z) { return M3(V3(X, 0.0f, 0.0f),
                                                   V3(0.0f, Y, 0.0f),
                                                   V3(0.0f, 0.0f, Z)); }
inline m3 M3Scale(v3 Dim) { return M3Scale(Dim.x, Dim.y, Dim.z); }
inline m4 M4Scale(f32 X, f32 Y, f32 Z, f32 W) { return M4(V4(X, 0.0f, 0.0f, 0.0f),
                                                          V4(0.0f, Y, 0.0f, 0.0f),
                                                          V4(0.0f, 0.0f, Z, 0.0f),
                                                          V4(0.0f, 0.0f, 0.0f, W)); }
inline m4 M4Scale(v3 Dim) { return M4Scale(Dim.x, Dim.y, Dim.z, 1.0f); }
inline m4 M4Scale(v4 Dim) { return M4Scale(Dim.x, Dim.y, Dim.z, Dim.w); }

// NOTE: Pos Matrix
inline m3 M3Pos(v2 Pos) { return M3(V3(1, 0, 0),
                                    V3(0, 1, 0),
                                    V3(Pos, 1)); }
inline m4 M4Pos(v3 Pos) { return M4(V4(1, 0, 0, 0),
                                    V4(0, 1, 0, 0),
                                    V4(0, 0, 1, 0),
                                    V4(Pos, 1)); }

// NOTE: Translate Matrix
inline m3 M3Translate(m3 Mat, v2 Pos) { return M3(Mat.v[0], Mat.v[1], V3(Mat.v[2].xy + Pos, Mat.v[2].z)); }
inline m4 M4Translate(m4 Mat, v3 Pos) { return M4(Mat.v[0], Mat.v[1], Mat.v[2], V4(Mat.v[3].xyz + Pos, Mat.v[3].w)); }

// NOTE: Rotation Matrix
inline m3 M3Rotation(f32 Angle) { f32 CAngle = Cos(Angle); f32 SAngle = Sin(Angle); return M3(V3(CAngle, SAngle, 0.0f), V3(-SAngle, CAngle, 0.0f), V3(0.0f, 0.0f, 1.0f)); }
inline m4 M4Rotation(f32 AngleX, f32 AngleY, f32 AngleZ)
{
    m4 RotX = M4(V4(1.0f, 0.0f, 0.0f, 0.0f),
                 V4(0.0f, Cos(AngleX), Sin(AngleX), 0.0f),
                 V4(0.0f, -Sin(AngleX), Cos(AngleX), 0.0f),
                 V4(0.0f, 0.0f, 0.0f, 1.0f));
    m4 RotY = M4(V4(Cos(AngleY), 0.0f, -Sin(AngleY), 0.0f),
                 V4(0.0f, 1.0f, 0.0f, 0.0f),
                 V4(Sin(AngleY), 0.0f, Cos(AngleY), 0.0f),
                 V4(0.0f, 0.0f, 0.0f, 1.0f));
    m4 RotZ = M4(V4(Cos(AngleZ), Sin(AngleZ), 0.0f, 0.0f),
                 V4(-Sin(AngleZ), Cos(AngleZ), 0.0f, 0.0f),
                 V4(0.0f, 0.0f, 1.0f, 0.0f),
                 V4(0.0f, 0.0f, 0.0f, 1.0f));
    return RotX*RotY*RotZ;
}
inline m4 M4Rotation(v3 Rotate) { return M4Rotation(Rotate.x, Rotate.y, Rotate.z); }

// TODO: These should be renamed since theyre mixed up, also call them flip not invert
inline m3 InvertXAxis()
{
    m3 Result = M3Identity();
    Result.e[0] = -1.0f;

    return Result;
}

inline m3 InvertYAxis()
{
    m3 Result = M3Identity();
    Result.v[1].y = -1.0f;

    return Result;
}

// =======================================================================================================================================
// NOTE: Common Aabb Functions 
// =======================================================================================================================================

inline aabb2 Enlarge(aabb2 Aabb, v2 AddRadius) { return AabbMinMax(Aabb.Min - AddRadius, Aabb.Max + AddRadius); }
inline aabb2i Enlarge(aabb2i Aabb, v2i AddRadius) { return AabbiMinMax(Aabb.Min - AddRadius, Aabb.Max + AddRadius); }
inline aabb3 Enlarge(aabb3 Aabb, v3 AddRadius) { return AabbMinMax(Aabb.Min - AddRadius, Aabb.Max + AddRadius); }

inline aabb2 Translate(aabb2 Aabb, v2 Displacement) { return AabbMinMax(Aabb.Min + Displacement, Aabb.Max + Displacement); }
inline aabb2i Translate(aabb2i Aabb, v2i Displacement) { return AabbiMinMax(Aabb.Min + Displacement, Aabb.Max + Displacement); }
inline aabb3 Translate(aabb3 Aabb, v3 Displacement) { return AabbMinMax(Aabb.Min + Displacement, Aabb.Max + Displacement); }

inline v2 AabbGetCenter(aabb2 A) { return Lerp(A.Min, A.Max, 0.5f); }
inline v2i AabbGetCenter(aabb2i A) { return Lerp(A.Min, A.Max, 0.5f); }
inline v3 AabbGetCenter(aabb3 A) { return Lerp(A.Min, A.Max, 0.5f); }

inline v2 AabbGetDim(aabb2 A) { return A.Max - A.Min; }
inline v2i AabbGetDim(aabb2i A) { return A.Max - A.Min; }
inline v3 AabbGetDim(aabb3 A) { return A.Max - A.Min; }

inline v2 AabbGetRadius(aabb2 A) { return (A.Max - A.Min) * 0.5f; }
inline v2i AabbGetRadius(aabb2i A) { return (A.Max - A.Min) / 2; }
inline v3 AabbGetRadius(aabb3 A) { return (A.Max - A.Min) * 0.5f; }

inline v3 AabbGetBotPos(aabb3 A) { return V3(AabbGetCenter(A).xy, A.Min.z); }

inline b32 Intersect(aabb2 A, v2 B) { return (A.Min.x <= B.x && A.Max.x >= B.x &&
                                              A.Min.y <= B.y && A.Max.y >= B.y); }
inline b32 Intersect(aabb3 A, v3 B) { return (A.Min.x <= B.x && A.Max.x >= B.x &&
                                              A.Min.y <= B.y && A.Max.y >= B.y &&
                                              A.Min.z <= B.z && A.Max.z >= B.z); }

inline b32 Intersect(aabb2 A, v2 Pos, f32 Radius)
{
    v2 Closest = V2(Clamp(Pos.x, A.Min.x, A.Max.x), Clamp(Pos.y, A.Min.y, A.Max.y));
    return LengthSquared(Pos - Closest) <= Radius;
}

inline b32 Intersect(aabb2 A, aabb2 B) { return (!((B.Min.x > A.Max.x) || (B.Max.x < A.Min.x) ||
                                                   (B.Min.y > A.Max.y) || (B.Max.y < A.Min.y))); }
inline b32 Intersect(aabb3 A, aabb3 B) { return (!((B.Min.x > A.Max.x) || (B.Max.x < A.Min.x) ||
                                                   (B.Min.y > A.Max.y) || (B.Max.y < A.Min.y) ||
                                                   (B.Min.z > A.Max.z) || (B.Max.z < A.Min.z))); }

inline b32 IntersectNotInclusive(aabb2 A, aabb2 B) { return (!((B.Min.x >= A.Max.x) || (B.Max.x <= A.Min.x) ||
                                                               (B.Min.y >= A.Max.y) || (B.Max.y <= A.Min.y))); }
inline b32 IntersectNotInclusive(aabb3 A, aabb3 B) { return (!((B.Min.x >= A.Max.x) || (B.Max.x <= A.Min.x) ||
                                                               (B.Min.y >= A.Max.y) || (B.Max.y <= A.Min.y) ||
                                                               (B.Min.z >= A.Max.z) || (B.Max.z <= A.Min.z))); }

inline v2 GetNearestPointOnAabbToPoint(aabb2 A, v2 Pos)
{
    v2 Result = V2(Clamp(Pos.x, A.Min.x, A.Max.x), Clamp(Pos.y, A.Min.y, A.Max.y));
    return Result;
}

inline v3 GetNearestPointOnAabbToPoint(aabb3 A, v3 Pos)
{
    v3 Result = V3(Clamp(Pos.x, A.Min.x, A.Max.x), Clamp(Pos.y, A.Min.y, A.Max.y), Clamp(Pos.z, A.Min.z, A.Max.z));
    return Result;
}

// TODO: Rename to point not circle?
inline f32 DistBetweenAabbPointSquared(aabb2 A, v2 Pos)
{
    v2 Closest = GetNearestPointOnAabbToPoint(A, Pos);
    f32 Result = LengthSquared(Closest - Pos);
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
inline aabb2i LoadAabb2i(aabb2i_soa Soa, u32 Index) { aabb2i Result = AabbiMinMax(LoadV2i(Soa.Min, Index), LoadV2i(Soa.Max, Index)); return Result; }
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

// NOTE: Gameplay functions
inline q4 Q4FacingDirection(f32 Angle)
{
    // TODO: Remove adjustment
    f32 AdjustedAngle = Angle + (1.0f / 2.0f)*Pi32;
    q4 Result = Q4AxisAngle(V3(0, 0, 1), AdjustedAngle);

    return Result;
}
