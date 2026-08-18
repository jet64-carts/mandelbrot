// carts/mandelbrot/config.h
//
// MANDELBROT — the single source of truth for every tunable in the cart. A number that shapes
// what the explorer does lives here, named once and grouped by the module that owns it; the
// structural constants of a module (a lane width, a matrix) live beside the code they shape.
// The relationships between these numbers are pinned by _Static_assert at the foot of the file.
#ifndef MANDELBROT_CONFIG_H
#define MANDELBROT_CONFIG_H

// --- Display ---------------------------------------------------------------------------------
// The explorer draws the machine's full default screen, one palette index per pixel.
#define MB_SCREEN_WIDTH 1280
#define MB_SCREEN_HEIGHT 720
#define MB_SCREEN_PIXELS (MB_SCREEN_WIDTH * MB_SCREEN_HEIGHT)

// --- The view (view.c) -----------------------------------------------------------------------
// The home view frames the whole set with a margin: centred a little left of the origin, and
// tall enough to hold the set's full height (it reaches +-1.12i) in a 16:9 screen.
#define MB_HOME_CENTRE_RE (-0.65)
#define MB_HOME_CENTRE_IM 0.0
#define MB_HOME_HEIGHT 2.5  // the imaginary extent of the screen at home

// How far out a player may zoom beyond home, in octaves (factors of two), and how far the
// centre may wander from the origin. Together they keep every visible c within
// MB_VIEW_MAX_RADIUS of the origin, which the kernel's bailout relies on (see below).
#define MB_ZOOM_OUT_LIMIT_OCTAVES 2.0
#define MB_CENTRE_LIMIT 3
#define MB_VIEW_MAX_RADIUS 16

// THE CROSSOVER. Direct double iteration resolves a pixel only while the pixel spacing is a
// comfortable number of ulps of |c|: below 2^-MB_DIRECT_PRECISION_BITS of the view's magnitude the
// grid's columns would be only 2^(52 - MB_DIRECT_PRECISION_BITS) ulps apart, and the orbits' own
// rounding would begin to show. Deeper, every pixel iterates by perturbation against a reference
// orbit in the bignum (reference.h); design.md §9.5 measures where the switch is seamless. The
// magnitude is floored at a quarter, the nearest the set's boundary comes to the origin.
#define MB_DIRECT_PRECISION_BITS 40
#define MB_DIRECT_MAGNITUDE_FLOOR 0.25

// THE DEPTH LIMIT: a pixel spacing of 2^-880 (about 1.2e-265, a magnification of about 10^262 from
// home). It is where the bignum's BIG_FRAC_LIMBS (30 limbs, 960 bits) run out: a reference serving
// it needs 880 bits and its 64 guard bits (reference.h), and 16 bits are kept in hand (bignum.h has
// the centre's own count). Every per-pixel quantity, held in doubles, is still well inside the
// double's range there (design.md §9.6): delta_c and |delta| reach down to the spacing, 2^142 above
// the smallest normal (delta^2 may underflow harmlessly: it is the term a BLA step drops, and a
// validity test that underflows only refuses a block); dz/dc and a BLA table's B, at most about
// 2^40 / spacing (measured on the deep suite), stay below 2^920, 2^104 under the largest double. The
// deepest minibrot of the deep suite, 10^-254 across, is framed at 2^-852. The view may not zoom
// deeper.
#define MB_DEPTH_LIMIT_BITS 880
#define MB_DEPTH_LIMIT_LOG2 (-(double)MB_DEPTH_LIMIT_BITS)

// --- Iteration limits (view.c) ---------------------------------------------------------------
// The escape radius, squared. Far larger than 2, for two reasons: the smooth iteration count is
// accurate to O(1/R^2), and an escaped orbit grows monotonically for ever once |z| exceeds
// max(2, |c|) — so a speculative block that ends inside the circle cannot have escaped within
// it. That second reason is why MB_VIEW_MAX_RADIUS must stay well inside the radius.
#define MB_BAILOUT_RADIUS 1024

// The automatic iteration limit: a base at the home view, rising by a fixed amount for every
// octave of zoom, clamped. The player can scale it by powers of two (camera.c).
#define MB_ITERATIONS_BASE 600
#define MB_ITERATIONS_PER_OCTAVE 160
#define MB_ITERATIONS_MIN 64
#define MB_ITERATIONS_MAX 262144  // a reference orbit's capacity (reference.h): a centre reference
                                  // must cover the limit, and WRAM holds two of 2^18
#define MB_ITERATION_SCALE_MAX_STEPS 6  // the player may halve or double this many times

// --- Scheduler (scheduler.c) -----------------------------------------------------------------
// The values in this section that sit behind #ifndef are the cart's; the guards exist so the
// harness (tests/mandelbrot/build.sh) can run the cart on fewer cores of a shared machine, with
// smaller quanta in its soak, and at other tile shapes in its ablation lane. The cart itself
// never overrides them.

// Worker threads in the pool; the runtime holds any request to its own maximum of eight.
#ifndef MB_WORKER_COUNT
#define MB_WORKER_COUNT 8
#endif

// A tile is the unit of scheduling: one job computes a band of one tile. The shapes measured were
// within noise of one another; this one keeps a batch inside the pool's job ring while giving
// eight workers enough units to balance (design.md §6.3). Both sides must be multiples of the
// coarsest refinement stride.
#ifndef MB_TILE_WIDTH
#define MB_TILE_WIDTH 128
#endif
#ifndef MB_TILE_HEIGHT
#define MB_TILE_HEIGHT 48
#endif

// The work each draw hands the pool, in plain lane-iterations (the kernel's own cost unit; a
// derivative lane-iteration counts at its measured price, kernel.h). It is deterministic by
// construction — planned from counted costs, never from a clock — so the same inputs give the same
// frames on every machine. The web player runs its jobs on the main thread, so its quanta are a
// single core's share, not yet measured in a browser.
//   BATCH    what each draw at rest launches, to run while the host presents and runs the next
//            _update (a draw on the move launches the motion quantum).
//            Measured at three workers from a standing start, 32 M converges the suite in about half
//            the draws of 16 M and 13-20% sooner; at eight workers through the agent door (which
//            runs ticks back to back, so every batch lands in its tick) the gate's scenarios keep a
//            median whole-tick p95 of 3.3 ms or less, with no late tick (design.md §15.4).
//   RESTART  the synchronous budget of a restart at rest (a jump the player will look at): the
//            coarsest level is paid for first, whatever it costs, and finer levels get what is left.
//   MOTION   the synchronous budget of a draw whose view moved or whose limit changed: the strips
//            reuse could not carry are computed before the frame is shown. A restart while the view
//            still moves gets this too — it is shown for a tick.
// Sized on the reference machine (an Apple M4, 4 performance and 6 efficiency cores, about 7-10 G
// lane-iterations a second across eight workers when shared; design.md §6.4).
#ifndef MB_QUANTUM_NATIVE
#define MB_QUANTUM_NATIVE 32000000LL
#define MB_QUANTUM_WEB 4000000LL
#define MB_RESTART_QUANTUM_NATIVE 24000000LL
#define MB_RESTART_QUANTUM_WEB 3000000LL
#define MB_MOTION_QUANTUM_NATIVE 16000000LL
#define MB_MOTION_QUANTUM_WEB 2000000LL
#endif
//   RIM      at depth, where the draw computes nothing else (design.md §9.9), the synchronous budget
//            of the coarsest level of the FAR previews a zoom has stretched beyond the old picture
//            (§6.3): sampled before the frame is shown, so the rim a zoom out uncovers shows the new
//            view coarse rather than the old one's edge smeared across it for a tick. A quarter of
//            the motion quantum: the rim a tick uncovers is a band, and its coarsest level one
//            sample in 256.
#ifndef MB_RIM_QUANTUM_NATIVE
#define MB_RIM_QUANTUM_NATIVE (MB_MOTION_QUANTUM_NATIVE / 4)
#define MB_RIM_QUANTUM_WEB (MB_MOTION_QUANTUM_WEB / 4)
#endif

// Reuse serves a view change only while at least this share of the columns and of the rows can be
// carried (a zoom of two octaves in one tick carries a quarter); below it, the view restarts.
#define MB_MIN_REUSE_SHARE 0.25

// The certificates' margins, as multiples of the theorem's least radius (s / sqrt 2 for a cell of
// side s). The interior fill is exact wherever the theorem holds, so it needs none beyond the
// rounding of its estimate; the exterior's buys interpolation accuracy (design.md §7.4).
// The BLA tables' validity tolerance: a block of steps is taken linearly while the terms it drops
// stay within this fraction of the terms it keeps (design.md §9.3). Chosen by measurement, the
// largest whose pictures stay within a small fraction of a colour step of plain perturbation's
// (tests/mandelbrot/deep_test.c --sweep; the sweep is recorded in design.md §15).
#ifndef MB_BLA_EPSILON
#define MB_BLA_EPSILON 0x1p-38
#endif

#ifndef MB_EXTERIOR_MARGIN
#define MB_EXTERIOR_MARGIN 2.0
#endif
#ifndef MB_INTERIOR_MARGIN
#define MB_INTERIOR_MARGIN 1.05
#endif

// The cost assumed for a sample before anything in its tile has been measured, in
// lane-iterations. Only the planner reads it; a wrong guess costs balance, never correctness.
#define MB_DEFAULT_SAMPLE_COST 256.0

// --- Colour (palette.c) ----------------------------------------------------------------------
// One turn of the colour cycle, in ticks (120 per second). The palette is regenerated at a
// fractional phase each tick; no pixel is rewritten.
#define MB_CYCLE_PERIOD_TICKS 2400

// The interior's Copper gradient, top and bottom of the screen: a deep night blue settling to
// black, drawn for free behind every palette-index-0 pixel.
#define MB_INTERIOR_TOP_R 10
#define MB_INTERIOR_TOP_G 14
#define MB_INTERIOR_TOP_B 34
#define MB_INTERIOR_BOTTOM_R 0
#define MB_INTERIOR_BOTTOM_G 0
#define MB_INTERIOR_BOTTOM_B 2

// --- Camera (camera.c) -----------------------------------------------------------------------
// Zoom eases toward its target: each tick closes this fraction of the remaining distance (in
// octaves), which gives a half-life of about four ticks. A remainder below the snap threshold
// is closed outright, so a settled view is exactly its target.
#define MB_ZOOM_EASE_PER_TICK 0.16
#define MB_ZOOM_SNAP_OCTAVES 0.001
#define MB_WHEEL_OCTAVES_PER_NOTCH 0.5
#define MB_KEY_ZOOM_OCTAVES_PER_SECOND 1.5
#define MB_TRIGGER_ZOOM_OCTAVES_PER_SECOND 2.0

// Keyboard and stick panning, in screen heights per second at full deflection, eased so a tap
// neither jerks nor overshoots.
#define MB_PAN_SCREENS_PER_SECOND 0.6
#define MB_PAN_EASE_PER_TICK 0.2
#define MB_PAN_SNAP_PIXELS 0.01  // a dying pan slower than this, per tick, stops outright

// --- HUD (hud.c) -----------------------------------------------------------------------------
#define MB_HUD_FONT 3  // the project's monogram face
#define MB_HUD_TEXT_SIZE 16
#define MB_HUD_MARGIN 14
#define MB_HUD_LINE_GAP 4

// --- Relationships ---------------------------------------------------------------------------
_Static_assert(MB_SCREEN_WIDTH % MB_TILE_WIDTH == 0, "tiles must tile the screen's width");
_Static_assert(MB_SCREEN_HEIGHT % MB_TILE_HEIGHT == 0, "tiles must tile the screen's height");
_Static_assert(MB_BAILOUT_RADIUS >= 2 * MB_VIEW_MAX_RADIUS,
               "an escaped orbit must grow monotonically everywhere a player can look");
_Static_assert(MB_CENTRE_LIMIT < MB_VIEW_MAX_RADIUS, "the centre limit bounds the view");
_Static_assert(MB_ITERATIONS_MIN <= MB_ITERATIONS_BASE && MB_ITERATIONS_BASE <= MB_ITERATIONS_MAX,
               "the base iteration limit must lie inside its clamp");
_Static_assert(MB_WORKER_COUNT >= 1, "the scheduler needs at least one worker");

#endif // MANDELBROT_CONFIG_H
