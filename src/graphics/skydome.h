#pragma once

/* HDR Sky dome */
class SkyDome {
    vector<f32> sampler;
    i32 w, h, n;

   public:
    SkyDome() = default;
    SkyDome(const char* file_path);

    float3 sample_dir(float3 dir) const {
        u32 u = w * atan2f(dir.z, dir.x) * INV2PI - 0.5f;
        u32 v = h * acosf(dir.y) * INVPI - 0.5f;
        u32 i = u + v * w;
        return 0.65f * float3(sampler[i * 3], sampler[i * 3 + 1], sampler[i * 3 + 2]);
    }
};
