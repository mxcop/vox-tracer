#pragma once

#ifdef DEV

/**
 * @brief Update the development control GUI.
 */
extern void devgui_control();

#else

/**
 * @brief DOES NOTHING IF "DEV" IS UNDEFINED.
 */
void devgui_control(){};

#endif
