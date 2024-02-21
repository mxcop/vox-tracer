#pragma once

#include "graphics/camera.h"
#include "graphics/primitives/voxel-volume.h"
#include "graphics/light.h"

class Renderer : public TheApp {
   public:
    void init();
    u32 trace(const Ray& ray) const;
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

    unique_ptr<VoxelVolume> volume;
    Camera camera;
    vector<LightSource> lights;
};
