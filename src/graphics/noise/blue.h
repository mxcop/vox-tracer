#pragma once

/* R2 irrationals */
constexpr f32 R2 = 1.22074408460575947536f;
constexpr f32 R2X = 1.0f / R2;
constexpr f32 R2Y = 1.0f / (R2 * R2);
constexpr f32 R2Z = 1.0f / (R2 * R2 * R2);
constexpr f32 R2_2D = 1.32471795724474602596f;
constexpr f32 R2X_2D = 1.0f / R2_2D;
constexpr f32 R2Y_2D = 1.0f / (R2_2D * R2_2D);

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

    float2 sample_2d(u32 x, u32 y) const {
        x = x % w, y = y % h;
        u32 i = (x + y * w) * n;
        return float2(sampler[i + 0], sampler[i + 1]);
    }
};
