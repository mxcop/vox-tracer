#pragma once

#ifdef DEV

class Renderer;

/* Development namespace (only accessible if "DEV" is defined) */
namespace dev {
extern bool hide_devgui;
extern Renderer* renderer;
}  // namespace dev

#endif
