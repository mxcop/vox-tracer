#include "blue.h"

#include <stb_image.h>

BlueNoise::BlueNoise() {
    f32* data = stbi_loadf("assets/blue_noise_512.png", &w, &h, &n, 0);
    sampler = vector<f32>(data, data + (w * h * n));
}
