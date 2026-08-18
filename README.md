# Mandelbrot

The Mandelbrot set on the [Jet64 fantasy console](https://73v3.itch.io), explored live at
120 Hz. Drag the plane with the mouse, spin the wheel to dive toward the pointer, and the
picture follows: a new view is on screen, complete, in the tick you ask for it, and sharpens
to every last pixel a moment later — at home, or at a magnification of 10²⁶⁰.

## Controls

| | Mouse and keyboard | Gamepad |
|---|---|---|
| Pan | drag with the left button; arrows or W A S D | left stick |
| Zoom | the wheel, about the pointer; E in, Q out | right trigger in, left trigger out |
| Return home | R | START |
| Fly to a bookmark | 1 to 9: home, the valleys and spirals, minibrots at 10⁻²⁹ and 10⁻¹³⁰ | |
| Tour the bookmarks | T, or leave it alone for a minute; any input stops it | B |
| Style | C: Relief or Silk | X |
| Gradient | V: Classic, Ember or Glacier | R3 |
| Colour bands | `[` thinner, `]` thicker | |
| Iterations | `-` fewer, `=` more | LB fewer, RB more |
| Colour cycle | SPACE pauses and resumes | A |
| Status line | H | SELECT |
| Help | `/` | Y |
| Screenshot | P | |

## How it is made

**The look** (`look.c`) shades the frame from the field in one pass over the pool: hue on a scale
measured from the view itself (the logarithm of μ above the view's least, across its middle 96%,
three turns of the gradient), so home and 10⁻²⁵⁴ are both the whole gradient; every filament
drawn as a line by its distance estimate, taken from the gradient of μ; and the field lit as a
surface from the upper left. **Bookmarks fly** (`flight.c`): out, across and in, each leg eased,
exact to the last bit of the target's centre at any depth, slowed wherever the picture falls
behind, and any input takes the camera back. Not built yet: van Wijk and Nuij's optimal flight path,
the Ink style, the Orchid and Jade gradients, per-bookmark styles, and an iteration limit driven by
the picture.

**The kernel** (`kernel.c`) iterates a stream of points in small groups of lanes, several
groups at once, so the compiler's vectoriser turns every step into SIMD arithmetic in plain
ISO C. Blocks of iterations run speculatively, with no escape test; a group that overran is
restored and replayed checked. Its answer is the plain scalar loop's, bit for bit, and it runs
at about three times the scalar loop's speed on one core. The main cardioid, the period-2 bulb
and exactly periodic orbits are recognised as interior without iterating to the limit.

**The scheduler** (`scheduler.c`) spends a fixed quantum of counted work each frame, never a
measured time, so the same inputs draw the same frames on every machine. A new view computes
its coarsest level at once and then refines coarse-to-fine across the thread pool, the
costliest work first, while the host presents the frame.

**Motion** (`reuse.c`) never starts again. As the view pans or zooms, every column and row of
the old picture that lies within half a pixel of the new one (a pixel, in a zoom) is carried
across, and only the strips it cannot carry are computed — before the frame is shown, so every
moving frame is whole. At depth, where a pixel costs thousands of times more, they are computed
beside it, and meanwhile every pixel shows the nearest value the picture holds: a neighbour's a
pixel away, or, where a zoom has outrun the old picture, a coarse sample of the new view — the
coarsest gaps filled first. A pan comes to rest on the pixels it already holds; a zoom settles and
refines in place, and no pixel is ever replaced by a worse one.

**Depth** (`bignum.c`, `reference.c`, `deep.c`). Below about 10¹² a double can no longer name
a pixel. From there the view's centre is held in a fixed-point bignum of up to 960 bits, one
reference orbit is iterated in it, and every pixel follows that orbit in doubles as a small
perturbation — rebasing onto the reference's start whenever it strays, and skipping the long
stretches where the perturbation is linear with a table of bivariate linear approximations — 5
to 45 times faster than stepping, the deeper the more. The next reference is found and built while the current one still
serves — the nucleus of the minibrot on screen, by Newton's method in perturbation form — so a
continuous dive never waits, and the switch from direct iteration is measured to be seamless.

**The certificates** (`certify.c`) decide whole squares of the picture without iterating any
pixel in them. Two theorems give a sample a disc that provably holds none of the set's boundary —
outside the set from the escape's derivative, inside it from the attracting cycle's — and where
the four corners of a square own discs that cover it, no filament, however fine, can hide
inside: the square is filled, its smooth count interpolated from the corners.

**The colour** (`palette.c`, `look.c`) is a smooth iteration count mapped onto a cyclic
gradient mixed in the perceptual OKLab space, dithered with Roberts' R2 sequence and interleaved
gradient noise, and cycled by rewriting the palette rather than a single pixel. The compute shades
nothing: the field holds values, and the look alone turns them into colour. The interior is the
Copper's true-colour gradient.

**The arithmetic** uses only `+ - * /` and `sqrt`, with `#pragma STDC FP_CONTRACT OFF` in
every unit, and its own `log2` and `exp2`, so the workstation, a desktop export and a browser
draw the same picture to the bit.

`docs/mandelbrot/design.md` in the Jet64 tree tells the whole story, with the measurements
behind every choice.

## Verifying it

The harness lives with the console rather than with the cart — from the root of the Jet64
tree:

```bash
bash tests/mandelbrot/build.sh
```

It checks every unit for warnings and for ISO C, the elementary
functions against libm, the bignum and the perturbation kernels against an independent ground
truth computed in exact arithmetic around seven minibrots down to 10⁻²⁵⁴, the crossover from
direct iteration, the kernel against a scalar oracle bit for bit over six views, the
floating-point contract across compile lines, the camera, every pixel of the renderer against
the oracle under three thread arrangements — every certified pixel on the right side of the set,
within a twentieth of a palette entry — the picture in motion frame by frame, at home and deep,
and the whole cart under the address, undefined behaviour and thread sanitisers. `--bench` adds
the throughput, convergence and ablation measurements.

Made with Jet64. MIT licence — see [LICENSE](LICENSE).
