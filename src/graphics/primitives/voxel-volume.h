#pragma once

#include "graphics/rays/hit.h"
#include "morton.h"

/* 16 voxels per unit of space */
constexpr f32 DEFAULT_VOXEL_SCALE = 16.0f;

/* Use morton ordering for voxels (slightly faster) */
#define USE_MORTON 0
/* Use the AVX2 gather instructions to load voxels (worse on my AMD CPU) */
#define USE_AVX2_GATHER 1

struct VoxelVolume {
    f32 scale = DEFAULT_VOXEL_SCALE;
    /* Bounding box */
    union {
        f128 bmin;
        struct {
            float3 bbmin;
            f32 _;
        };
    };
    union {
        f128 bmax;
        struct {
            float3 bbmax;
            f32 _;
        };
    };
    union {
        f128 vsize;
        struct {
            float3 voxel_size;
            f32 _;
        };
    };
    std::vector<u8> voxels;

    VoxelVolume() = default;
    /**
     * @param pos Position of the voxel volume in world space.
     * @param size Size of the voxel volume in voxels.
     * @param scale The world space scale of 1 voxel. (default 16 voxels per unit)
     */
    VoxelVolume(float3 pos, int3 size, f32 scale = DEFAULT_VOXEL_SCALE)
        : scale(scale), bbmin(pos), bbmax(pos + float3(size) / scale), voxel_size(size) {
#if USE_MORTON
        /* Morton only works for perfect cubes */
        u32 bytes = pow(max(max(size.x, size.y), size.z), 3);
#else
        u32 bytes = size.x * size.y * size.z;
#endif
        voxels = std::vector<u8>();
        voxels.resize(bytes);

        /* Populate the voxel data */
        f32 rx = 1.0f / 128.0f;
        for (u32 z = 0; z < size.z; ++z) {
            const f32 fz = (f32)z / 128.0f;
            for (u32 y = 0; y < size.y; ++y) {
                const f32 fy = (f32)y / 128.0f;
                f32 fx = 0;
                for (u32 x = 0; x < size.x; ++x, fx += rx) {
                    const f32 noise = noise3D(fx, fy, fz);
                    const f32 noise2 = noise3D(fx + 0.33f, fy + 0.66f, fz + 0.99f);
#if USE_MORTON
                    const u64 i = morton_encode(x, y, z);
#else
                    const u64 i = (z * size.x * size.y) + (y * size.x) + x;
#endif
                    voxels[i] = noise > 0.09f ? (noise2 * 0xFF) : 0x00;
                }
            }
        }
    }

    /**
     * @brief Intersect the voxel volume with a ray.
     */
    HitInfo intersect(const Ray& ray) const;
    bool is_occluded(const Ray& ray) const;

    /**
     * @brief Intersect the voxel volume with a packet of 4 rays.
     */
    PacketHitInfo intersect(const RayPacket& packet) const;

   private:
    /* Fast ray to aabb check, using SSE and FMA. */
    inline f32 ray_vs_aabb(const Ray& ray) const {
        /* Idea to use fmsub to save 1 instruction came from
         * <http://www.joshbarczak.com/blog/?p=787> */
        const f128 ord = _mm_mul_ps(ray.o, ray.rd);
        const f128 t1 = _mm_fmsub_ps(bmin, ray.rd, ord);
        const f128 t2 = _mm_fmsub_ps(bmax, ray.rd, ord);

        /* Find the near and far intersection point */
        const f128 vmax4 = _mm_max_ps(t1, t2), vmin4 = _mm_min_ps(t1, t2);
        const f32 tmax = _min(_min(vmax4.m128_f32[0], vmax4.m128_f32[1]), vmax4.m128_f32[2]);
        /* The last max(0) here handles being inside of the AABB */
        const f32 tmin =
            _max(_max(_max(vmin4.m128_f32[0], vmin4.m128_f32[1]), vmin4.m128_f32[2]), 0.0f);
        return tmin <= tmax ? tmin : BIG_F32;
    }

    /* Fast ray packet to aabb check, using SSE */
    void ray4_vs_aabb(const RayPacket& packet, PacketHitInfo& out) const {
        f128 tmin = _mm_setzero_ps();
        f128 tmax = _mm_set_ps1(BIG_F32);

        for (u32 a = 0; a < 3; ++a) {
            #if 0
            const f128 aabb_min = _mm_set_ps1(packet.signs[a] > 0.0f ? bbmin[a] : bbmax[a]);
            const f128 aabb_max = _mm_set_ps1(packet.signs[a] > 0.0f ? bbmax[a] : bbmin[a]);

            const f128 dmin = _mm_mul_ps(_mm_sub_ps(aabb_min, packet.ro[a]), packet.ird[a]);
            const f128 dmax = _mm_mul_ps(_mm_sub_ps(aabb_max, packet.ro[a]), packet.ird[a]);

            tmin = _mm_max_ps(tmin, dmin);
            tmax = _mm_min_ps(tmax, dmax);
            #else

            // TODO: optimize function check more...
            const f128 t1 =
                _mm_mul_ps(_mm_sub_ps(_mm_set_ps1(bbmin[a]), packet.ro[a]), packet.ird[a]);
            const f128 t2 =
                _mm_mul_ps(_mm_sub_ps(_mm_set_ps1(bbmax[a]), packet.ro[a]), packet.ird[a]);

            tmin = _mm_max_ps(tmin, _mm_min_ps(t1, t2));
            tmax = _mm_min_ps(tmax, _mm_max_ps(t1, t2));
            //tmin = _mm_min_ps(_mm_max_ps(t1, tmin), _mm_max_ps(t2, tmin));
            //tmax = _mm_max_ps(_mm_min_ps(t1, tmax), _mm_min_ps(t2, tmax));
            #endif

            //float t1 = (box->min[d] - ray->origin[d]) * ray->dir_inv[d];
            //float t2 = (box->max[d] - ray->origin[d]) * ray->dir_inv[d];

            //tmin = max(tmin, min(t1, t2));
            //tmax = min(tmax, max(t1, t2));
        }
        // tmin = _mm_max_ps(_mm_setzero_ps(), tmin);

        /* Use a mask to remove non-intersections (tmin > tmax) */
        const f128 mask = _mm_cmple_ps(tmin, tmax);
        out.depth = _mm_blendv_ps(BIG_F128, tmin, mask);
        out.exit_t = tmax;
    }

    /* Fetch a voxel from this volumes voxel data. */
    __forceinline u8 fetch_voxel(const int3& idx, const int3& size) const {
#if USE_MORTON
        const u64 i = morton_encode(idx.x, idx.y, idx.z);
#else
        size_t i = ((size_t)idx.z * size.x * size.y) + ((size_t)idx.y * size.x) + idx.x;
#endif
        return voxels[i];
    }

    /* Fetch 4 voxels using the indices inside an SSE register. */
    __forceinline i128 fetch_voxels(const i128 indices) const { 
#if USE_AVX2_GATHER
        return gather_voxels(indices);
#else
        return _mm_set_epi32(voxels[indices.m128i_u32[3]], voxels[indices.m128i_u32[2]],
                             voxels[indices.m128i_u32[1]], voxels[indices.m128i_u32[0]]);
#endif
    }

    /* Gather 4 voxels from this volumes voxel data using SIMD. */
    __forceinline i128 gather_voxels(const i128 indices) const {
        const i128 u8mask = _mm_set1_epi32(0x000000FF);
        return _mm_and_epi32(_mm_i32gather_epi32((int*)voxels.data(), indices, sizeof(u8)), u8mask);
    }
};
