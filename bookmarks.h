// carts/mandelbrot/bookmarks.h
//
// BOOKMARKS — places worth the trip, one key each (1 to 9), which the camera flies to (flight.h):
// the lead's curation, from home through the valleys and spirals to two deep minibrots.
//
// Beneath them, the DEEP SUITE: the seven minibrots of the harness's deep views, from 10^29 down
// to 10^254, each framed four of its own sizes across (two of them are bookmarks too).
#ifndef MANDELBROT_BOOKMARKS_H
#define MANDELBROT_BOOKMARKS_H

#include "view.h"

#define BOOKMARK_COUNT 9
#define SUITE_COUNT 7

/** The view of bookmark `index` (0 is home, on key 1); false if there is none. */
bool bookmark_view(int index, View* out);

/** Bookmark `index`'s name, or "" if there is none. */
const char* bookmark_name(int index);

/** The framed view of deep suite minibrot `index` (0 is the shallowest); false if there is none. */
bool suite_view(int index, View* out);

/** The period of deep suite minibrot `index`, or 0 if there is none. */
int suite_period(int index);

#endif // MANDELBROT_BOOKMARKS_H
