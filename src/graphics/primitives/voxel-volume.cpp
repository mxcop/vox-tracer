#include "voxel-volume.h"

constexpr int MAX_STEPS = 256;

// TODO: Multi-resolution traversal!
HitInfo VoxelVolume::intersect(const Ray& ray) const {
    HitInfo hit = {};

    /* Early return if bounding box was not hit */
    hit.depth = ray_vs_aabb(ray);
    if (hit.depth == BIG_F32) return hit;

    /* Calculate the size of the voxel volume in voxels */
    int3 size = floori(voxel_size);

    /* Calculate the starting voxel index */
    float3 ray_pos = ((ray.origin + ray.dir * (hit.depth + 0.001f)) - bbmin) * scale;
    int3 idx = floori(ray_pos);

    /* Clamp the starting voxel index inside the voxel volume */
    idx = clamp(idx, int3(0, 0, 0), size - 1);

    /* Which direction each axis will step in -1 or 1 */
    float3 step = ray.sign_dir;
    /* Indicates how far we must move (in units of t) to equal the width of a voxel */
    float3 delta = fabs(ray.r_dir);
    /* Determine t at which the ray crosses the first voxel boundary */
    float3 tmax = ((float3(idx) - ray_pos) + fmaxf(step, float3(0.0f))) * ray.r_dir;

    u32 i = 0;
    f32 t = 0.0f;
    float3 normal = {};
    for (; i < MAX_STEPS; ++i) {
        /* Fetch the current voxel */
        u8 voxel = fetch_voxel(idx, size);
        if (voxel != 0x00) {
            hit.depth += t / scale;
            hit.normal = -normal * step;
            hit.albedo = float3((f32)voxel / 256.0f);
            hit.steps = i;
            return hit;
        }

        /* Amanatides & Woo */
        /* <http://www.cse.yorku.ca/~amana/research/grid.pdf> */
        if (tmax.x < tmax.y) {
            if (tmax.x < tmax.z) {
                ray_pos.x += step.x;
                idx.x += step.x;
                if (idx.x < 0 || idx.x >= size.x) break;
                normal = float3(1, 0, 0);
                t = tmax.x;
                tmax.x += delta.x;
            } else {
                ray_pos.z += step.z;
                idx.z += step.z;
                if (idx.z < 0 || idx.z >= size.z) break;
                normal = float3(0, 0, 1);
                t = tmax.z;
                tmax.z += delta.z;
            }
        } else {
            if (tmax.y < tmax.z) {
                ray_pos.y += step.y;
                idx.y += step.y;
                if (idx.y < 0 || idx.y >= size.y) break;
                normal = float3(0, 1, 0);
                t = tmax.y;
                tmax.y += delta.y;
            } else {
                ray_pos.z += step.z;
                idx.z += step.z;
                if (idx.z < 0 || idx.z >= size.z) break;
                normal = float3(0, 0, 1);
                t = tmax.z;
                tmax.z += delta.z;
            }
        }
    }

    /* No hit occured! */
    hit.depth = BIG_F32;
    hit.steps = i;
    return hit;
}

bool VoxelVolume::is_occluded(const Ray& ray) const {
    HitInfo hit = {};

    /* Early return if bounding box was not hit */
    hit.depth = ray_vs_aabb(ray);
    if (hit.depth == BIG_F32) return false;

    /* Calculate the size of the voxel volume in voxels */
    const int3 size = floori(voxel_size);

    /* Calculate the starting voxel index */
    const float3 ray_pos = ((ray.origin + ray.dir * (hit.depth + 0.001f)) - bbmin) * scale;
    int3 idx = floori(ray_pos);

    /* Clamp the starting voxel index inside the voxel volume */
    idx = clamp(idx, int3(0, 0, 0), size - 1);

    /* Which direction each axis will step in -1 or 1 */
    const float3 step = ray.sign_dir;
    /* Indicates how far we must move (in units of t) to equal the width of a voxel */
    const float3 delta = fabs(ray.r_dir);
    /* Determine t at which the ray crosses the first voxel boundary */
    float3 tmax = ((float3(idx) - ray_pos) + fmaxf(step, float3(0.0f))) * ray.r_dir;

    for (u32 i = 0; i < MAX_STEPS; ++i) {
        /* Fetch the current voxel */
        u8 voxel = fetch_voxel(idx, size);
        if (voxel != 0x00) {
            return true;
        }

        /* Amanatides & Woo */
        /* <http://www.cse.yorku.ca/~amana/research/grid.pdf> */
        if (tmax.x < tmax.y) {
            if (tmax.x < tmax.z) {
                idx.x += step.x;
                if (idx.x < 0 || idx.x >= size.x) break;
                tmax.x += delta.x;
                //if (tmax.x >= 1.0f) return false;
            } else {
                idx.z += step.z;
                if (idx.z < 0 || idx.z >= size.z) break;
                tmax.z += delta.z;
                //if (tmax.z >= 1.0f) return false;
            }
        } else {
            if (tmax.y < tmax.z) {
                idx.y += step.y;
                if (idx.y < 0 || idx.y >= size.y) break;
                tmax.y += delta.y;
                //if (tmax.y >= 1.0f) return false;
            } else {
                idx.z += step.z;
                if (idx.z < 0 || idx.z >= size.z) break;
                tmax.z += delta.z;
                //if (tmax.z >= 1.0f) return false;
            }
        }
    }
    return false;
}

static __forceinline i128 normal_encode(const i128 x, const i128 y, const i128 z, const u32 sx,
                                        const u32 sy) {
    const i128 msx = _mm_set1_epi32((i32)sx);
    const i128 msy = _mm_set1_epi32((i32)sy);
    const i128 fc = _mm_mul_epi32(_mm_mul_epi32(z, msx), msy);
    const i128 sc = _mm_add_epi32(_mm_mul_epi32(y, msx), x);
    return _mm_add_epi32(fc, sc);
    // return (z * sx * sy) + (y * sx) + x;
}

/* Amanatides & Woo */
/* <http://www.cse.yorku.ca/~amana/research/grid.pdf> */
PacketHitInfo VoxelVolume::intersect(const RayPacket& packet) const {
    PacketHitInfo hit = {};

    /* Early return if no rays hit the bounding box */
    ray4_vs_aabb(packet, hit);
    const f128 miss_mask = _mm_cmpeq_ps(hit.depth, BIG_F128);
    if (_mm_movemask_ps(miss_mask) == 0b1111) return hit;

    /* Voxels per unit of space, and it's reciprocal */
    const f128 ms = _mm_set_ps1(scale), is = _mm_set_ps1(1.0f / scale);

    /* Determine the ray entry point into the volume */
    const f128 entry_t = _mm_add_ps(hit.depth, _mm_set_ps1(0.0000000001f));
    const f128 entry_x = _mm_fmadd_ps(packet.rd_x, entry_t, packet.ro_x);
    const f128 entry_y = _mm_fmadd_ps(packet.rd_y, entry_t, packet.ro_y);
    const f128 entry_z = _mm_fmadd_ps(packet.rd_z, entry_t, packet.ro_z);

    /* Ray position (floating) */
    const f128 rp_x = _mm_mul_ps(_mm_sub_ps(entry_x, _mm_set_ps1(bbmin.x)), ms);
    const f128 rp_y = _mm_mul_ps(_mm_sub_ps(entry_y, _mm_set_ps1(bbmin.y)), ms);
    const f128 rp_z = _mm_mul_ps(_mm_sub_ps(entry_z, _mm_set_ps1(bbmin.z)), ms);

    /* Ray step, for moving the ray in voxel space (-1 or 1) */
    const i128 step_x = _mm_set1_epi32(packet.signs.x);
    const i128 step_y = _mm_set1_epi32(packet.signs.y);
    const i128 step_z = _mm_set1_epi32(packet.signs.z);

    /* Distance (in units of t) equal to the size of a voxel */
    const f128 delta_x = _mm_abs_ps(packet.ird_x);
    const f128 delta_y = _mm_abs_ps(packet.ird_y);
    const f128 delta_z = _mm_abs_ps(packet.ird_z);

    /* Ray index, the 3D index of the current voxel */
    const i128 lwb = _mm_set1_epi32(0);
    i128q ri;
    ri.x = _mm_floorclamp_ps(rp_x, lwb, _mm_set1_epi32(voxel_size.x - 1));
    ri.y = _mm_floorclamp_ps(rp_y, lwb, _mm_set1_epi32(voxel_size.y - 1));
    ri.z = _mm_floorclamp_ps(rp_z, lwb, _mm_set1_epi32(voxel_size.z - 1));
    // ri.x = _mm_cvtps_epi32(_mm_floor_ps(rp_x));
    // ri.y = _mm_cvtps_epi32(_mm_floor_ps(rp_y));
    // ri.z = _mm_cvtps_epi32(_mm_floor_ps(rp_z));

    /* Next position of the ray (in voxel space) */
    const f128 next_x = _mm_sub_ps(_mm_cvtepi32_ps(ri.x), rp_x);
    const f128 next_y = _mm_sub_ps(_mm_cvtepi32_ps(ri.y), rp_y);
    const f128 next_z = _mm_sub_ps(_mm_cvtepi32_ps(ri.z), rp_z);

    /* Get the starting offset for tmax, if this isn't done a visual glitch will occur */
    const f128 offset_x = _mm_max_ps(_mm_cvtepi32_ps(step_x), _mm_setzero_ps());
    const f128 offset_y = _mm_max_ps(_mm_cvtepi32_ps(step_y), _mm_setzero_ps());
    const f128 offset_z = _mm_max_ps(_mm_cvtepi32_ps(step_z), _mm_setzero_ps());

    /* "t" at which the ray crosses the next voxel boundary */
    f128 tmax_x = _mm_mul_ps(_mm_add_ps(next_x, offset_x), packet.ird_x);
    f128 tmax_y = _mm_mul_ps(_mm_add_ps(next_y, offset_y), packet.ird_y);
    f128 tmax_z = _mm_mul_ps(_mm_add_ps(next_z, offset_z), packet.ird_z);

    /* "t" at which the ray exits the volume */
    // const f128 exit_t = _mm_sub_ps(hit.exit_t, _mm_set_ps1(0.00001f));

    /* Mask that tracks status of the rays (done = 0xFFFFFFFF) */
    i128 status_mask = (i128&)miss_mask;  // (i128&)_mm_cmpge_ps(hit.depth, BIG_F128);
    i128 exit_mask = (i128&)miss_mask;
    constexpr u32 RAYS = 4;

    /* Make sure the next position of the rays are safe! */
    // const f128 tnext = _mm_fmadd_ps(_mm_min_ps(_mm_min_ps(tmax_x, tmax_y), tmax_z), is, entry_t);
    /* Mark rays outside the volume as done */
    // status_mask = _mm_or_epi32(status_mask, (i128&)_mm_cmpge_ps(tnext, exit_t));
    // status_mask = _mm_or_epi32(status_mask, (i128&)_mm_cmplt_ps(tnext, entry_t));

    /* Early exit if all rays are done */
    // if (_mm_movemask_ps((f128&)status_mask) == 0b1111) {
    //     return hit;
    // }

    /* Traverse */
    u32 s = 0;
    f128 t = _mm_setzero_ps();
    for (; s < MAX_STEPS; ++s) {
        /* Mark rays outside the volume as done */
        // status_mask = _mm_or_epi32(status_mask, (i128&)_mm_cmpge_ps(t, exit_t));
        // status_mask = _mm_or_epi32(status_mask, (i128&)_mm_cmplt_ps(t, entry_t));

        // status_mask = _mm_or_epi32(status_mask, _mm_cmplt_epi32(ri.x, _mm_set1_epi32(0)));
        // status_mask = _mm_or_epi32(status_mask, _mm_cmplt_epi32(ri.y, _mm_set1_epi32(0)));
        // status_mask = _mm_or_epi32(status_mask, _mm_cmplt_epi32(ri.z, _mm_set1_epi32(0)));
        // status_mask =
        //     _mm_or_epi32(status_mask, _mm_cmpgt_epi32(ri.x, _mm_set1_epi32(voxel_size.x - 1)));
        // status_mask =
        //     _mm_or_epi32(status_mask, _mm_cmpgt_epi32(ri.y, _mm_set1_epi32(voxel_size.y - 1)));
        // status_mask =
        //     _mm_or_epi32(status_mask, _mm_cmpgt_epi32(ri.z, _mm_set1_epi32(voxel_size.z - 1)));

        // ri.x = _mm_max_epi32(ri.x, _mm_set1_epi32(0));
        // ri.y = _mm_max_epi32(ri.y, _mm_set1_epi32(0));
        // ri.z = _mm_max_epi32(ri.z, _mm_set1_epi32(0));
        // ri.x = _mm_min_epi32(ri.x, _mm_set1_epi32(voxel_size.x - 1));
        // ri.y = _mm_min_epi32(ri.y, _mm_set1_epi32(voxel_size.y - 1));
        // ri.z = _mm_min_epi32(ri.z, _mm_set1_epi32(voxel_size.z - 1));

        /* Get the voxel indices */
#if USE_MORTON
        const i128 indices = morton_encode(ri.x, ri.y, ri.z);
#else
        const i128 indices = normal_encode(ri.x, ri.y, ri.z, voxel_size.x, voxel_size.y);
#endif

        /* Mark any rays outside the volume as done */
        exit_mask =
            _mm_or_epi32(exit_mask, _mm_or_epi32(_mm_cmpgt_epi32(indices, _mm_set1_epi32(voxels.size() - 1)),
                         _mm_cmplt_epi32(indices, _mm_set1_epi32(0))));
        //t = _mm_blendv_ps(t, BIG_F128, (f128&)exit_mask);
        status_mask = _mm_or_epi32(status_mask, exit_mask);

        /* Early exit if all rays are done */
        if (_mm_movemask_ps((f128&)status_mask) == 0b1111) {
            // const f128 t = _mm_min_ps(_mm_min_ps(tmax_x, tmax_y), tmax_z);
            hit.depth = _mm_blendv_ps(_mm_fmadd_ps(t, is, entry_t), BIG_F128, (f128&)exit_mask);
            break;
        }

        /* Mark rays that hit a voxel as done */
        // const i128 safe_indices = _mm_clamp_epi32(indices, 0, voxels.size() - 1);
        const i128 safe_indices = _mm_andnot_epi32(status_mask, indices);
        const i128 vox = fetch_voxels(safe_indices);
        status_mask = _mm_or_epi32(status_mask, _mm_cmpgt_epi32(vox, _mm_setzero_si128()));

        /* Find the smallest axis for each ray */
        const f128 min_yz = _mm_min_ps(tmax_y, tmax_z);
        const f128 min_zx = _mm_min_ps(tmax_z, tmax_x);
        const f128 min_xy = _mm_min_ps(tmax_x, tmax_y);
        /* And make a mask */
        const f128 cmp_x = _mm_andnot_ps((f128&)status_mask, _mm_cmplt_ps(tmax_x, min_yz));
        const f128 cmp_y = _mm_andnot_ps((f128&)status_mask, _mm_cmplt_ps(tmax_y, min_zx));
        const f128 cmp_z = _mm_andnot_ps((f128&)status_mask, _mm_cmplt_ps(tmax_z, min_xy));

        /* Record the current "t" of each ray */
        // const f128 tmin_x = _mm_and_ps(cmp_x, tmax_x);
        // const f128 tmin_y = _mm_and_ps(cmp_y, tmax_y);
        // const f128 tmin_z = _mm_and_ps(cmp_z, tmax_z);
        // const f128 tmin = _mm_or_ps(_mm_or_ps(tmin_x, tmin_y), tmin_z);
        const f128 tmin = _mm_min_ps(_mm_min_ps(tmax_x, tmax_y), tmax_z);
        //t = tmin;
        t = _mm_blendv_ps(tmin, t, (f128&)status_mask);
        // t = _mm_fmadd_ps(tmin, is, entry_t);

        /* Step to the next voxel index */
        ri.x = _mm_add_epi32(ri.x, _mm_and_epi32((i128&)cmp_x, step_x));
        ri.y = _mm_add_epi32(ri.y, _mm_and_epi32((i128&)cmp_y, step_y));
        ri.z = _mm_add_epi32(ri.z, _mm_and_epi32((i128&)cmp_z, step_z));

        /* Update the ray tmax */
        tmax_x = _mm_add_ps(tmax_x, _mm_and_ps(cmp_x, delta_x));
        tmax_y = _mm_add_ps(tmax_y, _mm_and_ps(cmp_y, delta_y));
        tmax_z = _mm_add_ps(tmax_z, _mm_and_ps(cmp_z, delta_z));

        /* TODO: do this "mul" on the tmax and delta once at the start? */
        /* WARN: this isn't 100% correct because this is actually "next_t" */
        // t = _mm_fmadd_ps(_mm_min_ps(_mm_min_ps(tmax_x, tmax_y), tmax_z), is, entry_t);
    }

    /* DEBUG: Save the step count */
    hit.steps = s;
    return hit;
}
