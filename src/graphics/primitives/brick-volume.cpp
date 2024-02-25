#include "brick-volume.h"

constexpr u32 MAX_STEPS = 256;

BrickVolume::BrickVolume(const float3 pos, const int3 size, const f32 vpu)
    : bbmin(pos), bbmax(pos + float3(size) / vpu), vsize(size), bsize(size >> 3), vpu(vpu) {
    /* Volume size MUST BE a multiple of 8 */
    assert(size.x % 8 == 0 && size.y % 8 == 0 && size.z % 8 == 0);

    /* Setup the voxels and bricks */
    voxels.resize((u64)vsize.x * vsize.y * vsize.z);
    brickmap.resize((u64)bsize.x * bsize.y * bsize.z);

    /* TEMP: Random bricks */
    u32 seed = 18942663u;
    u32 seed2 = 74928367u;
    u32 i = 0;
    for (Brick512& brick : brickmap) {
        if (RandomFloat(seed) < 0.2f) {
            /* Allocate voxel packets */
            brick.packets = make_unique_for_overwrite<u8[]>(64);
            //brick.brick_index = i;

            /* Randomly assign voxels */
            for (u16 i = 0; i < 64; i++) {
                brick.packets[i] = 0u;
                for (u8 b = 0; b < 8; b++) {
                    if (RandomFloat(seed2) < 0.5f) {
                        brick.packets[i] |= (0b1 << b);
                        brick.popcnt++;
                    }
                }
            }
        }
        i++;
    }

    //for (u32 z = 0; z < bsize.z; z++) {
    //    for (u32 y = 0; y < bsize.y; y++) {
    //        for (u32 x = 0; x < bsize.x; x++) {
    //            const float3 noise_pos = float3(x, y, z) * 8.0f / 128.0f;
    //            const u32 i = (z * bsize.x * bsize.y) + (y * bsize.x) + x;

    //            const f32 noise = noise3D(noise_pos.x, noise_pos.y, noise_pos.z);
    //            if (noise > 0.09f) {
    //                /* Allocate voxel packets */
    //                brickmap[i].packets = make_unique_for_overwrite<u8[]>(64);
    //                // brick.brick_index = i;

    //                /* Randomly assign voxels */
    //                for (u16 j = 0; j < 64; j++) {
    //                    brickmap[i].packets[j] = 0u;
    //                    for (u8 b = 0; b < 8; b++) {
    //                        if (RandomFloat(seed2) < /*0.05f*/1.0f) {
    //                            brickmap[i].packets[j] |= (0b1 << b);
    //                            brickmap[i].popcnt++;
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //    }
    //}
}

/**
 * @brief Intersect the volume bounding box. (tmin > tmax, means no intersection)
 */
void BrickVolume::intersect_bb(const Ray& ray, f32& tmin_out, f32& tmax_out) const {
    /* Idea to use fmsub to save 1 instruction came from
     * <http://www.joshbarczak.com/blog/?p=787> */
    const f128 ord = _mm_mul_ps(ray.o, ray.rd);
    const f128 t1 = _mm_fmsub_ps(bmin, ray.rd, ord);
    const f128 t2 = _mm_fmsub_ps(bmax, ray.rd, ord);

    /* Find the near and far intersection point */
    const f128 vmax4 = _mm_max_ps(t1, t2), vmin4 = _mm_min_ps(t1, t2);

    /* Set the 4th element to 0 */
    const f128 tmin4 = (f128&)_mm_slli_si128((i128&)vmin4, 4);

    /* Get the horizontal minimum and maximum "t" */
    tmax_out = _mm_hmin3_ps(vmax4), tmin_out = _mm_hmax_ps(tmin4);
}

/**
 * @brief Intersect the volume with a ray.
 */
HitInfo BrickVolume::intersect(const Ray& ray) const {
    HitInfo hit;

    /* Exit if the ray misses the volume */
    f32 tmin, tmax;
    intersect_bb(ray, tmin, tmax);
    if (tmin > tmax - 0.01f) {
        hit.depth = BIG_F32;
        return hit;
    }
    hit.depth = tmin;

    /* Volume entry position */
    const float3 entry = ((ray.origin + ray.dir * (hit.depth + 0.01f)) - bbmin) * (vpu * 0.125f);
    /* Clamp the entry point inside the volume grid */
    int3 cell = clamp(floori(entry), int3(0), bsize - 1);

    /* Which direction each axis will step in -1 or 1 */
    const float3 step = ray.sign_dir;
    /* Indicates how far we must move (in units of t) to equal the width of a voxel */
    const float3 delta = fabs(ray.r_dir);
    /* Determine t at which the ray crosses the first voxel boundary */
    float3 side = ((float3(cell) - entry) + fmaxf(step, float3(0))) * ray.r_dir;

    f32 t = 0.0f;
    for (hit.steps = 0; hit.steps < MAX_STEPS; ++hit.steps) {
        /* Fetch the active cell */
        const Brick512& brick = get_brick(cell);
        if (brick.popcnt > 0) {
            const f32 brick_entry_t = hit.depth + t / (vpu * 0.125f);
            if (traverse_brick(brick, cell, ray, brick_entry_t, hit)) {
                hit.depth += brick_entry_t;
                return hit;
            }
        }

        /* Amanatides & Woo */
        /* <http://www.cse.yorku.ca/~amana/research/grid.pdf> */
        if (side.x < side.y) {
            if (side.x < side.z) {
                cell.x += step.x;
                if (cell.x < 0 || cell.x >= bsize.x) break;
                hit.normal = float3(1, 0, 0);
                t = side.x;
                side.x += delta.x;
            } else {
                cell.z += step.z;
                if (cell.z < 0 || cell.z >= bsize.z) break;
                hit.normal = float3(0, 0, 1);
                t = side.z;
                side.z += delta.z;
            }
        } else {
            if (side.y < side.z) {
                cell.y += step.y;
                if (cell.y < 0 || cell.y >= bsize.y) break;
                hit.normal = float3(0, 1, 0);
                t = side.y;
                side.y += delta.y;
            } else {
                cell.z += step.z;
                if (cell.z < 0 || cell.z >= bsize.z) break;
                hit.normal = float3(0, 0, 1);
                t = side.z;
                side.z += delta.z;
            }
        }
    }

    /* No hit occured! */
    hit.depth = BIG_F32;
    return hit;
}

bool BrickVolume::traverse_brick(const Brick512& brick, const int3& pos, const Ray& ray,
                                const f32 entry_t, HitInfo& hit) const {
    /* Brick minimum position */
    const float3 bmin = bbmin + float3(pos) / (vpu * 0.125f);
    /* Brick entry position */
    const float3 entry = ((ray.origin + ray.dir * (entry_t)) - bmin) * vpu;
    /* Clamp the entry point inside the volume grid */
    int3 cell = clamp(floori(entry), 0, 8 - 1);

    /* Which direction each axis will step in -1 or 1 */
    const float3 step = ray.sign_dir;
    /* Indicates how far we must move (in units of t) to equal the width of a voxel */
    const float3 delta = fabs(ray.r_dir);
    /* Determine t at which the ray crosses the first voxel boundary */
    float3 side = ((float3(cell) - entry) + fmaxf(step, float3(0))) * ray.r_dir;

    f32 t = 0.0f;
    for (hit.steps; hit.steps < MAX_STEPS; ++hit.steps) {
        /* Fetch the active cell */
        const u8 voxel = get_voxel(brick, cell);
        if (voxel) {
            hit.depth += t / vpu;
            return true;
        }

        /* Amanatides & Woo */
        /* <http://www.cse.yorku.ca/~amana/research/grid.pdf> */
        if (side.x < side.y) {
            if (side.x < side.z) {
                cell.x += step.x;
                if (cell.x < 0 || cell.x >= 8) break;
                hit.normal = float3(1, 0, 0);
                t = side.x;
                side.x += delta.x;
            } else {
                cell.z += step.z;
                if (cell.z < 0 || cell.z >= 8) break;
                hit.normal = float3(0, 0, 1);
                t = side.z;
                side.z += delta.z;
            }
        } else {
            if (side.y < side.z) {
                cell.y += step.y;
                if (cell.y < 0 || cell.y >= 8) break;
                hit.normal = float3(0, 1, 0);
                t = side.y;
                side.y += delta.y;
            } else {
                cell.z += step.z;
                if (cell.z < 0 || cell.z >= 8) break;
                hit.normal = float3(0, 0, 1);
                t = side.z;
                side.z += delta.z;
            }
        }
    }

    /* No hit occured! */
    return false;
}
