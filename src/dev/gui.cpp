#include "gui.h"

#ifdef DEV

/* Internal development variables */
namespace inter {
float3 light_color = float3(1);
f32 radius = 0.1f;
f32 power = 16.0f;
}  // namespace inter

/**
 * @brief Update the development control GUI.
 */
void devgui_control() {
    if (dev::hide_devgui) return;

    if (ImGui::Begin("Control")) {
        /* Create sphere lights from the GUI */
        if (ImGui::CollapsingHeader("Sphere Light")) {
            ImGui::ColorEdit3("Color", &inter::light_color.x);
            ImGui::InputFloat("Radius", &inter::radius);
            ImGui::InputFloat("Power", &inter::power);
            if (ImGui::Button("Spawn on camera")) {
                const float3 origin = dev::renderer->camera.pos;
                dev::renderer->area_lights.emplace_back(origin, inter::radius, inter::light_color,
                                                        inter::power);
                dev::renderer->reset_accu();
            }
        }
    }
    ImGui::End();
}

#endif
