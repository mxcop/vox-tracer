void Renderer::init() {
    /* Try load the camera settings */
    FILE* f = fopen("camera.bin", "rb");
    if (f) {
        fread(&camera, 1, sizeof(Camera), f);
        fclose(f);
    }

    /* Create a voxel volume */
    volume = make_unique<VoxelVolume>(float3(0.0f, 0.0f, 0.0f), int3(128, 128, 128));
}

u32 Renderer::trace(const Ray& ray) const {
    HitInfo hit = volume->intersect(ray);

    /* Skybox color if the ray missed */
    if (hit.depth == BIG_F32) return 0xFF101010;

    /* Shoot a shadow ray */
    const float3 light_dir = normalize(float3(0.5f, 0.2f, 0.8f));
    const Ray shadow_ray = Ray(ray.origin + ray.dir * (hit.depth - 0.001f), light_dir);
    bool in_shadow = volume->is_occluded(shadow_ray);

    float4 en = float4(hit.albedo, 1.0f);
    if (in_shadow) {
        en *= 0.2f;
    }
    return RGBF32_to_RGB8(&en);

    /* Visualize the distance */
    // u32 cd = fminf(hit.depth / 32.0f, 1.0f) * 0xFF;
    /* Visualize the step count */
    // u32 cd = (hit.steps / 256.0f) * 0xFF;
    // u32 color = (cd << 0) | (cd << 8) | (cd << 16) | 0xFF000000;
    // float4 en = float4(hit.albedo * dot(hit.normal, normalize(float3(0.5f, 0.5f, 0.5f))), 1.0f);
    // u32 color = RGBF32_to_RGB8(&en);
    // return color;
}

void Renderer::tick(f32 dt) {
    Timer t;

#if 0
    constexpr i32 TILE_SIZE = 16;
    /* Split the window into tiles */
#pragma omp parallel for schedule(dynamic, 1)
    for (i32 t = 0; t < (WIN_WIDTH * WIN_HEIGHT) / (TILE_SIZE * TILE_SIZE); ++t) {
        u32 sx = (t / (WIN_HEIGHT / TILE_SIZE)) * TILE_SIZE;
        u32 sy = (t % (WIN_HEIGHT / TILE_SIZE)) * TILE_SIZE;
        for (u32 y = 0; y < TILE_SIZE; ++y) {
            for (u32 x = 0; x < TILE_SIZE; ++x) {
                Ray ray = camera.get_primary_ray(sx + x, sy + y);
                u32 color = trace(ray);
                screen->pixels[(sx + x) + (sy + y) * WIN_WIDTH] = color;
            }
        }
    }
#elif 1
#pragma omp parallel for schedule(dynamic)
    for (i32 y = 0; y < WIN_HEIGHT; ++y) {
        for (i32 x = 0; x < WIN_WIDTH; ++x) {
            Ray ray = camera.get_primary_ray(x, y);
            u32 color = trace(ray);
            screen->pixels[x + y * WIN_WIDTH] = color;
        }
    }
#elif 0
    constexpr u32 TILE_SIZE = 16; // 7.7M rays/s
    for (u32 y = 0; y < WIN_HEIGHT; y += TILE_SIZE) {
        for (u32 x = 0; x < WIN_WIDTH; x += TILE_SIZE) {
            for (u32 v = 0; v < TILE_SIZE; ++v) {
                u32 yv = y + v;
                for (u32 u = 0; u < TILE_SIZE; ++u) {
                    u32 xu = x + u;
                    Ray ray = camera.get_primary_ray(xu, yv);
                    u32 color = trace(ray);
                    screen->pixels[xu + yv * WIN_WIDTH] = color;
                }
            }
        }
    }
#else
    // 12.6M rays/s
//#pragma omp parallel for schedule(dynamic)
    for (i32 y = 0; y < WIN_HEIGHT; y += 2) {
        for (i32 x = 0; x < WIN_WIDTH; x += 2) {
            const RayPacket packet = camera.get_primary_packet(x, y);
            const PacketHitInfo hit = volume->intersect(packet);

            for (u32 v = 0; v < 2; ++v) {
                for (u32 u = 0; u < 2; ++u) {
                    const f32 depth = hit.depth.m128_f32[u + v * 2];
                    /*if (depth == 0.0f) {
                        screen->pixels[(x + u) + (y + v) * WIN_WIDTH] = 0xFFFF1010;
                        continue;
                    }*/
                    if (depth >= BIG_F32) {
                        screen->pixels[(x + u) + (y + v) * WIN_WIDTH] = 0xFF101010;
                        continue;
                    }
                    const u32 cd = fminf(depth / 32.0f, 1.0f) * 0xFF;
                    // const u32 cd = (hit.steps / 256.0f) * 0xFF;
                    const u32 color = (cd << 0) | (cd << 8) | (cd << 16) | 0xFF000000;
                    screen->pixels[(x + u) + (y + v) * WIN_WIDTH] = color;
                }
            }
        }
    }
#endif

    frame_time = t.elapsed();

    /* Update the camera */
    camera.update(dt);
}

void Renderer::gui(f32 dt) {
    /* Window position */
    constexpr float PADDING = 10.0f;
    ImVec2 work_pos = ImGui::GetMainViewport()->WorkPos;
    ImGui::SetNextWindowPos(ImVec2(work_pos.x + PADDING, work_pos.y + PADDING), ImGuiCond_Always,
                            ImVec2(0.0f, 0.0f));

    constexpr ImGuiWindowFlags overlay_flags =
        ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_AlwaysAutoResize |
        ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing |
        ImGuiWindowFlags_NoNav | ImGuiWindowFlags_NoMove;

    /* Display performance */
    ImGui::SetNextWindowBgAlpha(0.35f);
    if (ImGui::Begin("Debug overlay", nullptr, overlay_flags)) {
        f32 frame_time_us = frame_time * 1'000'000.0f;
        f64 win_size = (f64)WIN_WIDTH * WIN_HEIGHT;
        ImGui::Text("Debug overlay\n");
        ImGui::Separator();
        ImGui::Text("FPS: %.1f", 1.0f / dt);
        ImGui::Separator();
        ImGui::Text("Ray/s: %.2fM", (win_size / frame_time) / 1'000'000.0);
        ImGui::Text("Ray time (mean): %.2fns", (frame_time_us / win_size) * 1'000.0);
        ImGui::Text("Ray time (goal): %.2fns", (0.0166666 / win_size) * 1.0e+9);
        ImGui::Separator();
        ImGui::Text("Frame time: %.2fms", frame_time * 1'000.0f);
    }
    ImGui::End();
}

void Renderer::shutdown() {
    /* Save the camera state */
    FILE* f = fopen("camera.bin", "wb");
    fwrite(&camera, 1, sizeof(Camera), f);
    fclose(f);
}