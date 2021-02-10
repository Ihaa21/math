
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
    Result = Transpose(Result);
    Result = Result*M4Pos(-Pos);
    
    return Result;
}

inline m4 LookAtM4(v3 View, v3 Up, v3 Right, v3 Pos)
{
    // IMPORTANT: We assume the provided vectors form a standard basis
    Assert(Abs(Dot(View, Right)) <= 0.0001f);
    Assert(Abs(Dot(Right, Up)) <= 0.0001f);
    Assert(Abs(Dot(View, Up)) <= 0.0001f);
    Assert(Abs(Length(View) - 1.0f) <= 0.0001f);
    Assert(Abs(Length(Up) - 1.0f) <= 0.0001f);
    Assert(Abs(Length(Right) - 1.0f) <= 0.0001f);

    m4 Result = M4Identity();
    Result.v[0].xyz = Right;
    Result.v[1].xyz = Up;
    Result.v[2].xyz = View;
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
    Result = Transpose(Result);
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
