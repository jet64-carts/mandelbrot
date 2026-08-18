// carts/mandelbrot/hud.c
#pragma STDC FP_CONTRACT OFF

#include "hud.h"

#include "config.h"
#include "numerics.h"
#include "palette.h"

#define HUD_LINE_BYTES 160
#define HUD_PERCENT 100.0

// Below this many octaves (a magnification of about 10 000), the magnification reads better as a
// plain multiple; beyond it, as a power of ten: "3.2 x 10^152".
#define HUD_PLAIN_OCTAVES 13.0
#define HUD_LOG10_2 0.30102999566398119521   // decimal digits an octave
#define HUD_LOG2_10 3.32192809488736234787   // octaves a decimal digit
#define HUD_DECIMAL 10.0
// A mantissa printed to one decimal: at or above this it would read "10.0", so it carries.
#define HUD_MANTISSA_CARRY 9.95

// The status line is outlined, not drop-shadowed: over the set's filaments a shadow on one side
// leaves the other edge of every stroke against a bright background. The outline is the text
// printed at these four offsets in the shadow ink, then once in its own.
static const int OUTLINE_OFFSETS[][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
#define OUTLINE_COUNT ((int)(sizeof(OUTLINE_OFFSETS) / sizeof(OUTLINE_OFFSETS[0])))

// The help card: a title, a table of controls in two columns of key and action, and a note for
// the gamepad beneath.
typedef struct {
	const char* key;
	const char* action;
} HelpEntry;

static const char HELP_TITLE[] = "Mandelbrot";
static const HelpEntry HELP_ENTRIES[] = {
	{"Drag", "pan"},
	{"Wheel", "zoom at the pointer"},
	{"Arrows, W A S D", "pan"},
	{"E / Q", "zoom in / out"},
	{"R", "return home"},
	{"1 - 9", "fly to a bookmark"},
	{"T", "tour the bookmarks"},
	{"- / =", "fewer / more iterations"},
	{"C", "change the style"},
	{"V", "change the gradient"},
	{"[ / ]", "thinner / thicker colour bands"},
	{"Space", "pause the colour"},
	{"H", "hide the status line"},
	{"P", "save a screenshot"},
	{"/", "close this card"},
};
#define HELP_ENTRY_COUNT ((int)(sizeof(HELP_ENTRIES) / sizeof(HELP_ENTRIES[0])))
#define HELP_COLUMNS 2
#define HELP_ROWS ((HELP_ENTRY_COUNT + HELP_COLUMNS - 1) / HELP_COLUMNS)
static const char* const HELP_PAD_LINES[] = {
	"Gamepad: the left stick pans, the triggers zoom, START returns home, LB and RB set",
	"the iterations, X changes the style, R3 the gradient, B tours, A pauses the colour, SELECT hides the HUD.",
};
#define HELP_PAD_LINE_COUNT ((int)(sizeof(HELP_PAD_LINES) / sizeof(HELP_PAD_LINES[0])))
#define HELP_PADDING 24
#define HELP_KEY_GAP 16     // between a key and its action
#define HELP_COLUMN_GAP 40  // between the two columns
#define HELP_BORDER 1
#define HELP_SECTION_LINES 1  // a blank line between title, table and note

static TextFormat text_format(uint8_t ink) {
	return jet_make_text_format(MB_HUD_FONT, ink, MB_HUD_TEXT_SIZE, 0, 0, 1.0f, 0, 0);
}

static void print_outlined(const char* text, int x, int y, uint8_t ink) {
	const TextFormat outline = text_format(PAL_HUD_SHADOW);
	for (int i = 0; i < OUTLINE_COUNT; i++) {
		jet_print_ext(text, x + OUTLINE_OFFSETS[i][0], y + OUTLINE_OFFSETS[i][1], &outline);
	}
	const TextFormat face = text_format(ink);
	jet_print_ext(text, x, y, &face);
}

// The magnification 2^octaves as mantissa x 10^exponent, the mantissa in [1, 10) once rounded to
// the one decimal it is printed with. octaves >= HUD_PLAIN_OCTAVES > 0.
static void power_of_ten(double octaves, double* mantissa, int* exponent) {
	const double decimal = octaves * HUD_LOG10_2;
	int e = (int)decimal;  // positive, so the conversion is the floor
	double m = num_exp2((decimal - (double)e) * HUD_LOG2_10);
	if (m >= HUD_MANTISSA_CARRY) {
		m /= HUD_DECIMAL;
		e++;
	}
	*mantissa = m;
	*exponent = e;
}

static void draw_status(const HudStatus* status) {
	char zoom[HUD_LINE_BYTES];
	if (status->octaves < HUD_PLAIN_OCTAVES) {
		snprintf(zoom, sizeof zoom, "%.1fx", num_exp2(status->octaves));
	} else {
		double mantissa;
		int exponent;
		power_of_ten(status->octaves, &mantissa, &exponent);
		snprintf(zoom, sizeof zoom, "%.1f x 10^%d", mantissa, exponent);
	}
	char line[HUD_LINE_BYTES];
	snprintf(line, sizeof line, "Zoom %s    Iterations %lld    Computed %.0f%%%s", zoom,
	         (long long)status->max_iterations, status->progress * HUD_PERCENT,
	         status->at_depth_limit ? "    Depth limit" : "");
	const TextFormat face = text_format(PAL_STATUS);
	const int y = MB_SCREEN_HEIGHT - MB_HUD_MARGIN - jet_text_line_height(&face);
	print_outlined(line, MB_HUD_MARGIN, y, PAL_STATUS);
}

// The widest key and the widest action of one column of the table.
static void column_widths(const TextFormat* body, int column, int* key_w, int* action_w) {
	*key_w = *action_w = 0;
	for (int i = column; i < HELP_ENTRY_COUNT; i += HELP_COLUMNS) {
		const int k = jet_text_width(HELP_ENTRIES[i].key, body);
		const int a = jet_text_width(HELP_ENTRIES[i].action, body);
		if (k > *key_w) *key_w = k;
		if (a > *action_w) *action_w = a;
	}
}

static void draw_help(void) {
	const TextFormat title = text_format(PAL_HUD_ACCENT);
	const TextFormat keys = text_format(PAL_HUD_ACCENT);
	const TextFormat body = text_format(PAL_HUD_TEXT);
	const TextFormat note = text_format(PAL_HUD_DIM);
	const int line_height = jet_text_line_height(&body) + MB_HUD_LINE_GAP;

	int key_x[HELP_COLUMNS], action_x[HELP_COLUMNS];
	int table_w = 0;
	for (int col = 0; col < HELP_COLUMNS; col++) {
		int key_w, action_w;
		column_widths(&body, col, &key_w, &action_w);
		key_x[col] = table_w;
		action_x[col] = table_w + key_w + HELP_KEY_GAP;
		table_w = action_x[col] + action_w + (col + 1 < HELP_COLUMNS ? HELP_COLUMN_GAP : 0);
	}
	int content_w = table_w;
	for (int i = 0; i < HELP_PAD_LINE_COUNT; i++) {
		const int w = jet_text_width(HELP_PAD_LINES[i], &note);
		if (w > content_w) content_w = w;
	}
	const int lines = 1 + HELP_SECTION_LINES + HELP_ROWS + HELP_SECTION_LINES + HELP_PAD_LINE_COUNT;
	const int card_w = content_w + 2 * HELP_PADDING;
	const int card_h = lines * line_height + 2 * HELP_PADDING;
	const int x = (MB_SCREEN_WIDTH - card_w) / 2;
	const int y = (MB_SCREEN_HEIGHT - card_h) / 2;
	jet_rectfill(x - HELP_BORDER, y - HELP_BORDER, card_w + 2 * HELP_BORDER, card_h + 2 * HELP_BORDER,
	             PAL_HUD_DIM);
	jet_rectfill(x, y, card_w, card_h, PAL_HUD_PANEL);

	const int left = x + HELP_PADDING;
	int row_y = y + HELP_PADDING;
	jet_print_ext(HELP_TITLE, left, row_y, &title);
	row_y += (1 + HELP_SECTION_LINES) * line_height;
	for (int i = 0; i < HELP_ENTRY_COUNT; i++) {
		const int col = i % HELP_COLUMNS;
		const int ey = row_y + (i / HELP_COLUMNS) * line_height;
		jet_print_ext(HELP_ENTRIES[i].key, left + key_x[col], ey, &keys);
		jet_print_ext(HELP_ENTRIES[i].action, left + action_x[col], ey, &body);
	}
	row_y += (HELP_ROWS + HELP_SECTION_LINES) * line_height;
	for (int i = 0; i < HELP_PAD_LINE_COUNT; i++) {
		jet_print_ext(HELP_PAD_LINES[i], left, row_y + i * line_height, &note);
	}
}

// The title card: the name, large, and one line of how to begin, centred a little above the middle.
#define TITLE_NAME_SIZE 48
#define TITLE_LINE_GAP 12
#define TITLE_HEIGHT_SHARE 0.4  // the card's top, as a share of the screen's height
static const char TITLE_NAME[] = "Mandelbrot";
static const char TITLE_LINE[] = "Drag to pan, wheel to zoom, 1 to 9 to travel, / for help";

static void draw_title(void) {
	const TextFormat name = jet_make_text_format(MB_HUD_FONT, PAL_TITLE, TITLE_NAME_SIZE, 0, 0, 1.0f, 0, 0);
	const TextFormat line = text_format(PAL_TITLE);
	const int top = (int)(MB_SCREEN_HEIGHT * TITLE_HEIGHT_SHARE);
	const int name_w = jet_text_width(TITLE_NAME, &name);
	const int line_w = jet_text_width(TITLE_LINE, &line);
	const TextFormat name_outline = jet_make_text_format(MB_HUD_FONT, PAL_HUD_SHADOW, TITLE_NAME_SIZE, 0, 0, 1.0f, 0, 0);
	for (int i = 0; i < OUTLINE_COUNT; i++) {
		jet_print_ext(TITLE_NAME, (MB_SCREEN_WIDTH - name_w) / 2 + OUTLINE_OFFSETS[i][0], top + OUTLINE_OFFSETS[i][1],
		              &name_outline);
	}
	jet_print_ext(TITLE_NAME, (MB_SCREEN_WIDTH - name_w) / 2, top, &name);
	print_outlined(TITLE_LINE, (MB_SCREEN_WIDTH - line_w) / 2, top + jet_text_line_height(&name) + TITLE_LINE_GAP,
	               PAL_TITLE);
}

void hud_draw(const HudStatus* status, bool show_status, bool show_help, bool show_title) {
	if (show_title && !show_help) draw_title();
	if (show_status) draw_status(status);
	if (show_help) draw_help();
}
