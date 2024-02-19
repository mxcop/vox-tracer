#pragma once

/* Contains 4, 3 component vectors, stored in vertical fashion. */
union f128q {
    f128 q[3] = {};
    struct {
        f128 x, y, z;
    };
};

/* Contains 4, 3 component vectors, stored in vertical fashion. */
union i128q {
    i128 q[3] = {};
    struct {
        i128 x, y, z;
    };
};

/* SSE, Clamp and then floor. */
inline static i128 _mm_floorclamp_ps(const f128 m, const i128 min, const i128 max) {
    return _mm_min_epi32(_mm_max_epi32(_mm_cvtps_epi32(_mm_floor_ps(m)), min), max);
}

__forceinline static i128 _mm_clamp_epi32(const i128 v, const u32 min, const u32 max) {
    return _mm_min_epi32(_mm_max_epi32(v, _mm_set1_epi32(min)), _mm_set1_epi32(max));
}

/* SSE, Absolute float vector. */
inline static f128 _mm_abs_ps(const f128 m) { return _mm_andnot_ps(_mm_set_ps1(-0.0f), m); }
