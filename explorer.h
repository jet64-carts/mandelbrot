// carts/mandelbrot/explorer.h
//
// THE EXPLORER — the four hooks' one conductor.
//
// Each tick is model, compute, shade and present, kept apart:
//   _update   controls -> camera (the view) and the colour cycle's phase: all state that keeps
//             time lives here, advanced once per tick;
//   _draw     fence the previous batch (and restart on a new view) -> shade the field into the
//             look's frame and present it with the palette at its phase and the HUD over it ->
//             launch the next batch.
// Rendering compute runs only from _draw, so a player that runs several _updates per presented
// frame never runs several batches. _init computes the home picture's restart, so the first
// frame is already refined.
#ifndef MANDELBROT_EXPLORER_H
#define MANDELBROT_EXPLORER_H

#include "jet64.h"

/**
 * Sets the screen, reserves every buffer, starts the workers and computes the home picture's
 * restart. False if WRAM ran out.
 */
bool explorer_init(void);

/** One tick of the model. */
void explorer_update(void);

/** One presented frame. */
void explorer_draw(void);

/** Stops the workers. */
void explorer_end(void);

#endif // MANDELBROT_EXPLORER_H
