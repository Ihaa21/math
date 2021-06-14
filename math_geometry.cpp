
// =======================================================================================================================================
// NOTE: aabb2
// =======================================================================================================================================

// NOTE: Aabb functions

inline aabb2 AabbMinMax(v2 Min, v2 Max) { return { Min, Max }; }
inline aabb2 AabbCenterRadius(v2 Center, v2 Dim) { return { Center - Dim, Center + Dim }; }

// =======================================================================================================================================
// NOTE: aabb2i
// =======================================================================================================================================

inline aabb2i Aabb2i(aabb2 AabbF) { return { V2i(AabbF.Min), V2i(AabbF.Max) }; }
inline aabb2i AabbMinMax(v2i Min, v2i Max) { return { Min, Max }; }
inline aabb2i AabbCenterRadius(v2i Center, v2i Dim) { return { Center - Dim, Center + Dim }; }

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
// NOTE: Lines (Infinite)
// =======================================================================================================================================

inline b32 LinesIntersect(v2 Start0, v2 Dir0, v2 Start1, v2 Dir1, v2* OutIntersection)
{
    b32 Result = false;
    f32 Epsilon = 1.0f / 1024.0f;
    
    // NOTE: Convert to slope intercept form
    f32 Slope1 = Abs(Dir0.x) < Epsilon ? NAN : (Dir0.y / Dir0.x);
    f32 Slope2 = Abs(Dir1.x) < Epsilon ? NAN : (Dir1.y / Dir1.x);
    
    // NOTE: Check for parallel lines
    if (!(Abs(Slope1 - Slope2) < Epsilon))
    {
        Result = true;
        if (IsNan(Slope1) && !IsNan(Slope2))
        {
            OutIntersection->x = Start0.x;
            OutIntersection->y = (Start0.x - Start1.x) * Slope2 + Start1.y;
        }
        else if (!IsNan(Slope1) && IsNan(Slope2))
        {
            OutIntersection->x = Start1.x;
            OutIntersection->y = (Start1.x - Start0.x) * Slope1 + Start1.y;
        }
        else
        {
            OutIntersection->x = (Slope1 * Start0.x - Slope2 * Start1.x + Start1.y - Start0.y) / (Slope1 - Slope2);
            OutIntersection->y = Slope2 * (OutIntersection->x - Start1.x) + Start1.y;
        }
    }
    
    return Result;
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
// NOTE: Common Aabb Functions 
// =======================================================================================================================================

inline aabb2 Enlarge(aabb2 Aabb, v2 AddRadius) { return AabbMinMax(Aabb.Min - AddRadius, Aabb.Max + AddRadius); }
inline aabb2i Enlarge(aabb2i Aabb, v2i AddRadius) { return AabbMinMax(Aabb.Min - AddRadius, Aabb.Max + AddRadius); }
inline aabb3 Enlarge(aabb3 Aabb, v3 AddRadius) { return AabbMinMax(Aabb.Min - AddRadius, Aabb.Max + AddRadius); }

inline aabb2 Translate(aabb2 Aabb, v2 Displacement) { return AabbMinMax(Aabb.Min + Displacement, Aabb.Max + Displacement); }
inline aabb2i Translate(aabb2i Aabb, v2i Displacement) { return AabbMinMax(Aabb.Min + Displacement, Aabb.Max + Displacement); }
inline aabb3 Translate(aabb3 Aabb, v3 Displacement) { return AabbMinMax(Aabb.Min + Displacement, Aabb.Max + Displacement); }

inline v2 AabbGetCenter(aabb2 A) { return Lerp(A.Min, A.Max, V2(0.5f)); }
inline v2i AabbGetCenter(aabb2i A) { return Lerp(A.Min, A.Max, V2(0.5f)); }
inline v3 AabbGetCenter(aabb3 A) { return Lerp(A.Min, A.Max, V3(0.5f)); }

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
