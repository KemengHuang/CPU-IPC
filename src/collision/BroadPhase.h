#pragma once

// Compile-time switch between the legacy uniform-grid SpatialHash and the new
// BVH-based broad-phase. Define CIPC_USE_BVH (e.g. via CMake option) to select
// BVHBroadPhase; otherwise SpatialHash is used. Callers should use the
// BroadPhase alias so the same source compiles with either implementation.

#include "collision/Ground.h"

#ifdef CIPC_USE_BVH
#include "collision/BVHBroadPhase.h"
using BroadPhase = BVHBroadPhase;
#else
#include "collision/SpatialHash.h"
using BroadPhase = SpatialHash;
#endif
