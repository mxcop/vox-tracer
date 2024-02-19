#pragma once

/* Basic SIMD friendly Ray structure. */
struct Ray {
    /* Ray origin */
    union {
        f128 o;
        struct {
            float3 origin;
            f32 _;
        };
    };
    /* Ray direction */
    union {
        f128 d;
        struct {
            float3 dir;
            f32 _;
        };
    };
    /* Ray direction reciprocal */
    union {
        f128 rd;
        struct {
            float3 r_dir;
            f32 _;
        };
    };
    /* Ray direction sign (-1 | 1) */
    union {
        f128 sd;
        struct {
            float3 sign_dir;
            f32 _;
        };
    };

    Ray() = default;
    inline Ray(const float3& origin, const float3& dir)
        : origin(origin), dir(dir), r_dir(1.0f / dir), sign_dir(sign_of_dir(dir)) {}

   private:
    /* Get the signs of a vector (-1 | 1) */
    inline float3 sign_of_dir(float3 d) const {
        float3 signs = float3();
        signs.x = copysign(1.0f, dir.x);
        signs.y = copysign(1.0f, dir.y);
        signs.z = copysign(1.0f, dir.z);
        return signs;
    }
};
