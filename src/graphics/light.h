#pragma once

struct LightSource {
    /* The amount of light expressed with 3 color channels */
    float3 light = float3(1.0f, 1.0f, 1.0f);
    float3 origin;

    LightSource(float3 o, float3 l) : light(l), origin(o) {}
};
