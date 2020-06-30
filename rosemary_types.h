#if !defined(ROSEMARY_TYPES_H)
/* ========================================================================
   $File: $
   $Date: $
   $Revision: $
   $Creator: Ihor Szlachtycz $
   $Notice: (C) Copyright 2014 by Dream.Inc, Inc. All Rights Reserved. $
   ======================================================================== */

#include <stdint.h>
#include <stddef.h>
#include <float.h>
#include <atomic>

typedef int8_t i8;
typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;

typedef std::atomic<i8> atomic_i8;
typedef std::atomic<i16> atomic_i16;
typedef std::atomic<i32> atomic_i32;
typedef std::atomic<i64> atomic_i64;

#define I16_MIN -32768
#define I16_MAX 32767
#define I32_MIN INT_MIN
#define I32_MAX INT_MAX

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef std::atomic<u8> atomic_u8;
typedef std::atomic<u16> atomic_u16;
typedef std::atomic<u32> atomic_u32;
typedef std::atomic<u64> atomic_u64;

#define U32_MAX 0xFFFFFFFF

typedef float f32;
typedef double f64;

typedef std::atomic<f32> atomic_f32;
typedef std::atomic<f64> atomic_f64;

#define F32_MAX FLT_MAX
#define F32_MIN -F32_MAX
#define F32_EPSILON FLT_EPSILON

typedef size_t mm;
typedef uintptr_t uptr;
typedef intptr_t iptr;

typedef std::atomic<mm> atomic_mm;
typedef std::atomic<uptr> atomic_uptr;
typedef std::atomic<iptr> atomic_iptr;

typedef uint8_t b8;
typedef uint16_t b16;
typedef uint32_t b32;

typedef std::atomic<b8> atomic_b8;
typedef std::atomic<b16> atomic_b16;
typedef std::atomic<b32> atomic_b32;

#define internal static
#define global static
#define local_global static

#if __clang__
#define Assert(Expression) if (!(Expression)) {__builtin_trap();}
#else
#define Snprintf _snprintf_s
#define Assert(Expression) if (!(Expression)) {*(int*)0 = 0;}
#endif
#define InvalidCodePath Assert(!"Invalid Code Path")

#define ArrayCount(Array) (sizeof(Array) / sizeof((Array)[0]))
// TODO: Can we get the sizeof the dereferenced value?
//#define GetArrayIndex(Base, Elem) (((u64)Elem - (u64)Base) / sizeof(Elem))

#define KiloBytes(Val) ((Val)*1024LL)
#define MegaBytes(Val) (KiloBytes(Val)*1024LL)
#define GigaBytes(Val) (MegaBytes(Val)*1024LL)
#define TeraBytes(Val) (GigaBytes(Val)*1024LL)

// NOTE: This is here because visual studio sucks at expanding VA_ARGS https://stackoverflow.com/questions/5134523/msvc-doesnt-expand-va-args-correctly
#define EXPAND( x ) x

#define ROSEMARY_TYPES_H
#endif
