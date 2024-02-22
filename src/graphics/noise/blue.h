#pragma once

//#define STB_IMAGE_IMPLEMENTATION
//#include "stb_image.h"

/* Blue noise sampler */
class BlueNoise {
    vector<f32> sampler;
    i32 w, h, n;

   public:
    BlueNoise();

    float3 sample_3d(u32 x, u32 y) const {
        x = x % w, y = y % h;
        u32 i = (x + y * w) * n;
        return float3(sampler[i + 0], sampler[i + 1], sampler[i + 2]);
    }
};
