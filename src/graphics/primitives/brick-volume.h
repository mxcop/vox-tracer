#pragma once

/* 16 voxels per unit of world space */
constexpr f32 DEFAULT_VPU = 16.0f;

class BrickVolume {
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
    /* Size of the volume in voxels */
    const int3 vsize;
    /* Size of the volume in bricks */
    const int3 bsize;

    /* Voxels per unit of world space */
    const f32 vpu = DEFAULT_VPU;
    /* Bricks per unit of world space */
    const f32 bpu = DEFAULT_VPU * 0.125f;

    /* Voxel data */
    // vector<u8> voxels;

    /* 8x8x8 Brick contains 512 voxels */
    struct Brick512 {
        /* Pointer to 64 voxel packets each packet contains 8 voxels */
        u8* packets = nullptr;
        /* Number of active voxels in the brick */
        u16 popcnt = 0;
        /* TODO: Have 6 more bytes here to use... */

        Brick512() = default;
        ~Brick512() {
            if (packets) delete[] packets;
        }
        Brick512(const Brick512&) = delete;
        Brick512(Brick512&&) = default;
        Brick512& operator=(const Brick512&) = delete;
        Brick512& operator=(Brick512&&) = default;
    };
    /* Brick map */
    Brick512* brickmap = nullptr;

    /**
     * @brief Intersect the volume bounding box. (tmin > tmax, means no intersection)
     */
    __forceinline void intersect_bb(const Ray& ray, f32& tmin_out, f32& tmax_out) const;

    /* Get a brick from the brickmap. */
    __forceinline const Brick512* get_brick(const int3& i) const {
        const u32 idx = (i.z * bsize.x * bsize.y) + (i.y * bsize.x) + i.x;
        return brickmap + idx;
    };
    /* Get a voxel from a brick. */
    __forceinline u8 get_voxel(const Brick512* brick, const int3& i) const {
        // u64 idx = morton_encode(i.x, i.y, i.z);
        const u32 idx = (i.z * 8 * 8) + (i.y * 8) + i.x;
        const u32 byte = idx >> 3;
        const u8 bitmask = 0b1 << (idx - (byte << 3));
        return brick->packets[byte] & bitmask;
    };

    __inline f32 traverse_brick(const Brick512* brick, const int3& pos, const Ray& ray,
                                 const f32 entry_t, HitInfo& hit) const;
    __inline bool traverse_brick(const Brick512* brick, const int3& pos, const Ray& ray,
                                 const f32 entry_t) const;

   public:
    BrickVolume() = default;
    BrickVolume(const float3 pos, const int3 size, const f32 vpu = DEFAULT_VPU);
    ~BrickVolume();

    BrickVolume(const BrickVolume&) = delete;
    BrickVolume(BrickVolume&&) = default;
    BrickVolume& operator=(const BrickVolume&) = delete;
    BrickVolume& operator=(BrickVolume&&) = default;

    /**
     * @brief Intersect the volume with a ray.
     */
    HitInfo intersect(const Ray& ray) const;
    /**
     * @brief Intersect and return true if the ray hit anything.
     */
    bool is_occluded(const Ray& ray) const;
};
