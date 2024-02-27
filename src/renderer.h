#pragma once

#include "graphics/camera.h"
#include "graphics/primitives/voxel-volume.h"
#include "graphics/primitives/brick-volume.h"
#include "graphics/light.h"
#include "graphics/noise/blue.h"
#include "graphics/skydome.h"

class Renderer : public TheApp {
   public:
    void init();
    u32 trace(const Ray& ray, const u32 x, const u32 y) const;
    void tick(f32 dt);
    void gui(f32 dt);
    void shutdown();
    /* User input */
    void MouseUp(int button) {}
    void MouseDown(int button) {}
    void MouseMove(int x, int y) { mousePos.x = x, mousePos.y = y; }
    void MouseWheel(float y) {}
    void KeyUp(int key) {}
    void KeyDown(int key) {}

    int2 mousePos;
    f32 frame_time = 1.0f;
    u32 frame = 0u;
    float3 sun_dir = {-0.619501f, 0.465931f, -0.631765f};

    unique_ptr<VoxelVolume> volume;
    //unique_ptr<BrickVolume> volume;
    Camera camera;
    vector<LightSource> lights;

    BlueNoise bnoise;
    SkyDome skydome;

    /* Accumulator */
    float4* accu = nullptr;
    mutable u32 accu_len = 1u;
};
