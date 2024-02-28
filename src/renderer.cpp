#include <functional>
#include <graphics/tonemap.h>

void printb(u64 u) {
    for (int i = 63; i >= 0; i--) {
        if ((i + 1) % 8 == 0) printf("_");
        printf("%i", (u & ((u64)1u << i)) ? 1 : 0);
    }
    printf("\n");
}
void printb(u32 u) {
    for (int i = 0; i < 32; i++) {
        if (i % 8 == 0) printf("_");
        printf("%i", u & (1 << i) ? 1 : 0);
    }
    printf("\n");
}

void Renderer::init() {
    /* Try load the camera settings */
    FILE* f = fopen("camera.bin", "rb");
    if (f) {
        fread(&camera, 1, sizeof(Camera), f);
        fclose(f);
    }

    /* Create the accumulator */
    accu = (float4*)MALLOC64(WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    if (accu) memset(accu, 0, WIN_WIDTH * WIN_HEIGHT * sizeof(float4));

    bnoise = BlueNoise();
    skydome = SkyDome("assets/kiara_1_dawn_8k.hdr");

    /* Create a voxel volume */
    //volume = new VoxelVolume(float3(0.0f, 0.0f, 0.0f), int3(128, 128, 128));
    volume = new BrickVolume(float3(0.0f, 0.0f, 0.0f), int3(2048, 2048, 2048));
}

/* Source : <https://github.com/tqjxlm/Monte-Carlo-Ray-Tracer> */
static float3 sample_hemisphere_weighted(const f32 r1, const f32 r2, const float3& n) {
    f32 theta = acos(sqrt(1.0f - r1));
    f32 phi = 2.0f * PI * r2;
    f32 xs = sinf(theta) * cosf(phi);
    f32 ys = cosf(theta);
    f32 zs = sinf(theta) * sinf(phi);
    float3 h = n;

    if ((abs(h.x) <= abs(h.y)) && (abs(h.x) <= abs(h.z))) {
        h.x = 1.0;
    } else if ((abs(h.y) <= abs(h.x)) && (abs(h.y) <= abs(h.z))) {
        h.y = 1.0;
    } else {
        h.z = 1.0;
    }

    float3 x = normalize(cross(h, n));
    float3 z = normalize(cross(x, n));

    return normalize(xs * x + ys * n + zs * z);
}

float3 sample_hemisphere_uniform(const float3& n) {
    float3 R;
    do {
        R = float3(RandomFloat() * 2 - 1, RandomFloat() * 2 - 1, RandomFloat() * 2 - 1);
    } while (dot(R, R) > 1.0f);
    R = normalize(R);
    if (dot(n, R) < 0.0f) {
        return -R;
    } else {
        return R;
    }
}

u32 Renderer::trace(const Ray& ray, const u32 x, const u32 y) const {
    const HitInfo hit = volume->intersect(ray);

    float4 color = float4(0);

/* Point & spot lights */
#if 0
    /* Skybox color if the ray missed */
    if (hit.depth >= BIG_F32) {
        color = skydome.sample_dir(ray.dir);
        return RGBF32_to_RGB8(&color);
        // return 0xFF101010;
    }

    const float3 hit_pos = ray.origin + ray.dir * hit.depth + hit.normal * 0.0001f;

    /* R2 irrationals */
    const f32 R2 = 1.22074408460575947536f;
    const f32 R2X = 1.0f / R2;
    const f32 R2Y = 1.0f / (R2 * R2);
    const f32 R2Z = 1.0f / (R2 * R2 * R2);
    const f32 R2_2D = 1.32471795724474602596f;
    const f32 R2X_2D = 1.0f / R2_2D;
    const f32 R2Y_2D = 1.0f / (R2_2D * R2_2D);

    /* Ambient light */
#if 1
    
#define CONSINE_SAMPLING 1

    float3 ambient_c = float3(0.0f);
    constexpr f32 SAMPLES = 1;
    for (u32 i = 0; i < SAMPLES; i++) {
        /* Generate a random direction */
#if CONSINE_SAMPLING
        /* Blue noise + R2 (cosine weighted distribution) */
        const float2 raw_noise = bnoise.sample_2d(x, y);
        const f32 quasi_x = fmod(raw_noise.x + R2X_2D * (f32)(frame + i), 1.0f);
        const f32 quasi_y = fmod(raw_noise.y + R2Y_2D * (f32)(frame + i), 1.0f);
        const float3 ambient_dir = sample_hemisphere_weighted(quasi_x, quasi_y, hit.normal);
#else
        /* Blue noise + R2 (uniform distribution) */
        // float3 raw_noise = bnoise.sample_3d(x, y);
        // float3 quasi_noise;
        // quasi_noise.x = fmod(raw_noise.x + R2X * (f32)(frame + i), 1.0f);
        // quasi_noise.y = fmod(raw_noise.y + R2Y * (f32)(frame + i), 1.0f);
        // quasi_noise.z = fmod(raw_noise.z + R2Z * (f32)(frame + i), 1.0f);
        // float3 ambient_dir = normalize(quasi_noise);
        float3 ambient_dir = sample_hemisphere_uniform(hit.normal);
#endif

        /* Shoot the ambient ray */
        const Ray ambient_ray = Ray(hit_pos, ambient_dir * 32.0f);
        const bool in_shadow = volume->is_occluded(ambient_ray);

        if (not in_shadow) {
            /* Multiply by 2 Pi because the area we're integrating is 2 Pi */
            const float3 sample = skydome.sample_dir(ambient_dir);
            ambient_c += (sample * TWOPI);
        }
    }
    /* Divide by the number of samples */
    color += hit.albedo * (ambient_c * (1.0f / SAMPLES));
    // color += ambient_c * 0.05f; /* TODO: remove this */

#endif

    /* Directional light */
    for (u32 i = 0; i < 1; i++) {
        /* Jitter light position (soft shadows) */
        const f32 JITTER_DIAMETER = 6.0f / 16.0f;
        const f32 JITTER_RADIUS = JITTER_DIAMETER * 0.5f;

        float3 raw_noise = bnoise.sample_3d(x, y);
        u32 frame_offset = frame % 128u;
        float3 quasi_noise;
        quasi_noise.x = fmod(raw_noise.x + R2X * (f32)frame_offset, 1.0f);
        quasi_noise.y = fmod(raw_noise.y + R2Y * (f32)frame_offset, 1.0f);
        quasi_noise.z = fmod(raw_noise.z + R2Z * (f32)frame_offset, 1.0f);
        const float3 jitter = (quasi_noise * JITTER_DIAMETER - JITTER_RADIUS);
        const float3 sun_dirj = normalize(sun_dir + jitter);

        /* Do nothing if the normal faces away from the light */
        f32 incidence = dot(hit.normal, sun_dirj);
        if (incidence <= 0.0f) continue;

        /* Shoot shadow ray */
        const Ray shadow_ray = Ray(hit_pos, sun_dirj * 32.0f);
        const bool in_shadow = volume->is_occluded(shadow_ray);

        if (not in_shadow) {
            const float3 sun_light = float3(3.0f, 2.5f, 2.0f);
            color += hit.albedo * sun_light * incidence;
        }
    }

    for (const LightSource& light : lights) {
#if 1
        /* Jitter light position (soft shadows) */
        const f32 JITTER_DIAMETER = 3.0f / 16.0f;
        const f32 JITTER_RADIUS = JITTER_DIAMETER * 0.5f;

        float3 raw_noise = bnoise.sample_3d(x, y);
        u32 frame_offset = frame % 128u;
        float3 quasi_noise;
        quasi_noise.x = fmod(raw_noise.x + R2X * (f32)frame_offset, 1.0f);
        quasi_noise.y = fmod(raw_noise.y + R2Y * (f32)frame_offset, 1.0f);
        quasi_noise.z = fmod(raw_noise.z + R2Z * (f32)frame_offset, 1.0f);
        const float3 light_pos = light.origin + (quasi_noise * JITTER_DIAMETER - JITTER_RADIUS);
#else
        const float3 light_pos = light.origin;
#endif

        const f32 light_dist = length(light_pos - hit_pos);
        const float3 light_dir = (light_pos - hit_pos) / light_dist;

        /* Do nothing if the normal faces away from the light */
        f32 incidence = dot(hit.normal, light_dir);
        if (incidence <= 0.0f) continue;

        /* Handle aperture of spot light */
        if (light.aperture < 1.0f) {
            const f32 a = (dot(light_dir, light.dir) + 1.0f) * 0.5f;
            if (a > light.aperture) continue;
            incidence *= 1.0f - (a / light.aperture);
        }

        const Ray shadow_ray = Ray(light_pos, hit_pos - light_pos);
        const bool in_shadow = volume->is_occluded(shadow_ray);

        /* Do nothing if the point is in shadow */
        if (in_shadow) continue;
        const f32 sqd = light_dist * light_dist;

        /* Area contribution */
        const float3 area_c = light.light * JITTER_DIAMETER;
        color += hit.albedo * area_c * incidence / sqd;
    }
#else
     color = float4((hit.normal + 1.0f) * 0.5f, 1.0f);
    //color = float4(hit.depth / 8.0f, hit.depth / 8.0f, hit.depth / 8.0f, 1.0f);
     //color = float4(hit.albedo, 1.0f);
    // color = float4(hit.steps / 64.0f, hit.steps / 64.0f, hit.steps / 64.0f, 1.0f);
     return RGBF32_to_RGB8(&color);
#endif

    /* Update accumulator */
    accu[x + y * WIN_WIDTH] += aces_approx(color);
    color = accu[x + y * WIN_WIDTH] / (f32)accu_len;

    return RGBF32_to_RGB8(&color);
}

void Renderer::tick(f32 dt) {
    frame++;
    if (frame > 100) frame = 0;
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
#elif 0
#pragma omp parallel for schedule(dynamic)
    for (i32 y = 0; y < WIN_HEIGHT; ++y) {
        for (i32 x = 0; x < WIN_WIDTH; ++x) {
            Ray ray = camera.get_primary_ray(x, y);
            u32 color = trace(ray, x, y);
            screen->pixels[x + y * WIN_WIDTH] = color;
        }
    }
#elif 1
    constexpr u32 TILE_SIZE = 16;  // 7.7M rays/s
#pragma omp parallel for schedule(dynamic)
    for (i32 y = 0; y < WIN_HEIGHT; y += TILE_SIZE) {
        for (u32 x = 0; x < WIN_WIDTH; x += TILE_SIZE) {
            for (u32 v = 0; v < TILE_SIZE; ++v) {
                u32 yv = y + v;
                for (u32 u = 0; u < TILE_SIZE; ++u) {
                    u32 xu = x + u;
                    Ray ray = camera.get_primary_ray(xu, yv);
                    u32 color = trace(ray, xu, yv);
                    screen->pixels[xu + yv * WIN_WIDTH] = color;
                }
            }
        }
    }
#else
    // 12.6M rays/s
#pragma omp parallel for schedule(dynamic)
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

    accu_len++;
    frame_time = t.elapsed();

    /* Update the camera */
    if (camera.update(dt)) {
        accu_len = 1u;
        memset(accu, 0, WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    }
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
    if (ImGui::Begin("Perf overlay", nullptr, overlay_flags)) {
        f32 frame_time_us = frame_time * 1'000'000.0f;
        f64 win_size = (f64)WIN_WIDTH * WIN_HEIGHT;
        ImGui::Text("Perf overlay\n");
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

    static float3 light_color = float3(1);
    static f32 aperture = 1.0f;
    if (ImGui::Begin("Debug window")) {
        ImGui::Text("Debug options\n");
        ImGui::Separator();
        //ImGui::InputFloat3("Cam dir", camera.target.cell);
        //ImGui::DragFloat3("Sun", &sun_dir.x, 0.01f, -1.0f, 1.0f);
        ImGui::Separator();
        ImGui::SliderFloat("Aperture", &aperture, 0, 1);
        ImGui::Separator();
        ImGui::ColorPicker3("New light", &light_color.x);
    }
    ImGui::End();

    static bool q_down = false;
    if (IsKeyDown(GLFW_KEY_Q) && q_down == false) {
        lights.emplace_back(camera.pos, normalize(camera.target - camera.pos), aperture,
                            light_color * 8.0f);
        q_down = true;
    }
    if (!IsKeyDown(GLFW_KEY_Q) && q_down == true) {
        q_down = false;
    }
}

void Renderer::shutdown() {
    /* Save the camera state */
    FILE* f = fopen("camera.bin", "wb");
    fwrite(&camera, 1, sizeof(Camera), f);
    fclose(f);
}