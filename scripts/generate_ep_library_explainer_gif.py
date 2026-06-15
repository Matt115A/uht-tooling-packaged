"""EP Library Profiler — explainer GIF.

Scenes:
  1. Sequence  — long reads fall onto the reference plasmid map
  2. Analyse   — per-position mismatch rate chart builds up left-to-right
  3. Model     — λ value counts up; six example gene copies show mutation loads
  4. Report    — Poisson variant-distribution histogram
"""

from __future__ import annotations
import math
import random
from pathlib import Path
from typing import List, Tuple

from PIL import Image, ImageDraw, ImageFilter, ImageFont

WIDTH, HEIGHT = 1200, 675
FPS          = 24
SCENE_FRAMES = 48
FADE_FRAMES  = 12
HOLD_FRAMES  = 6

OUTPUT_PATH = Path("assets/animations/ep_library_profiler.gif")
POSTER_PATH = Path("assets/animations/ep_library_profiler_poster.png")

# ── Colours ───────────────────────────────────────────────────────────────────
BG    = (245, 245, 247)
CARD  = (255, 255, 255)
INK   = (29,  29,  31)
MID   = (110, 110, 115)
BLUE  = (0,   113, 227)
BORD  = (210, 210, 215)
TEAL  = (0,   196, 175)
GREEN = (40,  185, 115)
AMBER = (255, 180, 50)
CORAL = (255, 80,  65)
VIOLET= (155, 100, 240)
SLATE = (75,  115, 160)

# ── Reference map fractions ───────────────────────────────────────────────────
ROI_START = 0.22
ROI_END   = 0.68

# ── Pre-computed read marks ───────────────────────────────────────────────────
_rnd = random.Random(42)
NUM_READS = 6
READ_MARKS: List[List[float]] = []
for _ri in range(NUM_READS):
    _m: List[float] = []
    for _ in range(_rnd.randint(3, 6)):
        _m.append(_rnd.uniform(ROI_START + 0.02, ROI_END - 0.02))
    if _rnd.random() < 0.55:
        _m.append(_rnd.uniform(0.03, ROI_START - 0.02))
    if _rnd.random() < 0.40:
        _m.append(_rnd.uniform(ROI_END + 0.02, 0.95))
    READ_MARKS.append(_m)

# ── Pre-computed mismatch profile ─────────────────────────────────────────────
PROFILE_LEN    = 700
PROFILE_MAX    = 0.028
BACKGROUND_RATE = 0.0042
_pr = random.Random(7)
_PROFILE: List[float] = []
for _pi in range(PROFILE_LEN):
    _f = _pi / PROFILE_LEN
    if ROI_START <= _f < ROI_END:
        _v = 0.016 + 0.006 * math.sin(_pi * 0.25) + _pr.gauss(0, 0.003)
    else:
        _v = 0.004 + _pr.gauss(0, 0.0008)
    _PROFILE.append(max(0.0, _v))

# ── Model parameters ──────────────────────────────────────────────────────────
LAMBDA_AA = 2.0

_COPY_MUTATIONS: List[List[float]] = [
    [],
    [0.40],
    [0.18, 0.67],
    [0.31, 0.52, 0.79],
    [0.25],
    [0.12, 0.43, 0.62, 0.84],
]
_COPY_LABELS = ["WT (0 mut.)", "1 mut.", "2 mut.", "3 mut.", "1 mut.", "4 mut."]

# ── Poisson histogram ─────────────────────────────────────────────────────────
_LAM = LAMBDA_AA
_hv  = [math.exp(-_LAM) * (_LAM ** k) / math.factorial(k) for k in range(4)]
_hv.append(max(0.0, 1.0 - sum(_hv)))
HIST_VALS   = _hv
HIST_LABELS = ["0", "1", "2", "3", "4+"]
HIST_COLS   = [TEAL, BLUE, VIOLET, AMBER, CORAL]

_FONT_PATH = "/System/Library/Fonts/Supplemental/Arial.ttf"
_BOLD_PATH = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"


def _font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    for path in ([_BOLD_PATH, _FONT_PATH] if bold else [_FONT_PATH, _BOLD_PATH]):
        try:
            return ImageFont.truetype(path, size)
        except OSError:
            continue
    return ImageFont.load_default()


F_HUGE = _font(68, bold=True)
F_BIG  = _font(44, bold=True)
F_MED  = _font(20)
F_SM   = _font(14)
F_TINY = _font(11, bold=True)

# ── Easing ────────────────────────────────────────────────────────────────────

def _ease_io(t: float) -> float:
    t = max(0.0, min(1.0, t))
    return 0.5 - 0.5 * math.cos(math.pi * t)

def _ease_out(t: float) -> float:
    t = max(0.0, min(1.0, t))
    return 1 - (1 - t) ** 3

def _elastic_out(t: float) -> float:
    if t <= 0: return 0.0
    if t >= 1: return 1.0
    p, s = 0.38, 0.095
    return 2 ** (-9 * t) * math.sin((t - s) * 2 * math.pi / p) + 1

# ── Primitive helpers ─────────────────────────────────────────────────────────

def _base_img() -> Image.Image:
    return Image.new("RGBA", (WIDTH, HEIGHT), BG)

def _draw_card(img: Image.Image, x0: int, y0: int, x1: int, y1: int) -> None:
    sh = Image.new("RGBA", img.size, (0, 0, 0, 0))
    sd = ImageDraw.Draw(sh)
    sd.rounded_rectangle((x0+4, y0+8, x1+4, y1+8), radius=18, fill=(0,0,0,32))
    img.alpha_composite(sh.filter(ImageFilter.GaussianBlur(12)))
    ImageDraw.Draw(img).rounded_rectangle(
        (x0, y0, x1, y1), radius=18,
        fill=(*CARD, 255), outline=(*BORD, 255), width=1)

def _scene_text(draw: ImageDraw.ImageDraw,
                title: str, subtitle: str, alpha: int) -> None:
    draw.rectangle((70, 488, WIDTH-70, 490), fill=(*BLUE, alpha//3))
    draw.text((80, 502),  title,    font=F_HUGE, fill=(*INK, alpha))
    draw.text((82, 588),  subtitle, font=F_MED,  fill=(*MID, alpha))

def _textw(font: ImageFont.ImageFont, s: str) -> int:
    try:
        return int(font.getlength(s))
    except AttributeError:
        return len(s) * 8

def _draw_ref_bar(img: Image.Image, x: int, y: int, w: int, h: int,
                  alpha: int = 255) -> None:
    draw = ImageDraw.Draw(img)
    draw.rounded_rectangle((x, y, x+w, y+h), radius=h//2,
                            fill=(248,248,250,alpha), outline=(*BORD,alpha), width=1)
    inner = w - 8
    # Backbone left
    b1w = int(inner * ROI_START)
    draw.rounded_rectangle((x+4, y+4, x+4+b1w, y+h-4),
                            radius=(h-8)//2, fill=(*SLATE, alpha))
    # ROI / gene
    rx  = x + 4 + b1w + 2
    rw  = int(inner * (ROI_END - ROI_START))
    draw.rounded_rectangle((rx, y+4, rx+rw, y+h-4),
                            radius=(h-8)//2, fill=(*GREEN, alpha))
    # Gene label inside ROI
    lbl = "gene"
    tw  = _textw(F_TINY, lbl)
    draw.text((rx + rw//2 - tw//2, y + h//2 - 6), lbl,
              font=F_TINY, fill=(255,255,255,alpha))
    # Backbone right
    b2x = rx + rw + 2
    b2w = (x + w - 4) - b2x
    if b2w > 2:
        draw.rounded_rectangle((b2x, y+4, b2x+b2w, y+h-4),
                                radius=(h-8)//2, fill=(*SLATE, alpha))

def _draw_read(img: Image.Image, x: int, y: int, w: int, h: int,
               alpha: int, marks: List[float]) -> None:
    draw = ImageDraw.Draw(img)
    draw.rounded_rectangle((x, y, x+w, y+h), radius=h//2,
                            fill=(215,222,232,alpha), outline=(*BORD,alpha), width=1)
    for frac in marks:
        px = x + int(w * frac)
        py = y + h // 2
        r  = 4
        draw.ellipse((px-r, py-r, px+r, py+r), fill=(*CORAL, alpha))

def _profile_to_y(rate: float, chart_top: int, chart_bottom: int) -> int:
    return int(chart_bottom - (rate / PROFILE_MAX) * (chart_bottom - chart_top))


# ═══════════════════════════════════════════════════════════════════════════════
# SCENE 1 — SEQUENCE
# ═══════════════════════════════════════════════════════════════════════════════

def scene_sequence(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)

    REF_X, REF_Y, REF_W, REF_H = 110, 60, 980, 28
    READ_H, READ_SP = 22, 54

    ref_a = int(255 * min(1.0, t * 8))
    _draw_ref_bar(img, REF_X, REF_Y, REF_W, REF_H, ref_a)

    draw = ImageDraw.Draw(img)
    if ref_a > 120:
        roi_cx = REF_X + int(REF_W * (ROI_START + ROI_END) / 2)
        draw.text((roi_cx - 18, REF_Y - 20), "ROI", font=F_SM, fill=(*GREEN, ref_a))
        draw.text((REF_X + 4, REF_Y - 20), "plasmid reference", font=F_SM, fill=(*MID, ref_a))

    for idx in range(NUM_READS):
        delay   = idx * 0.09
        local_t = max(0.0, min(1.0, (t - delay) / max(1e-4, 1.0 - delay)))
        anim_t  = _elastic_out(local_t)
        y_final = REF_Y + REF_H + 14 + idx * READ_SP
        y       = int(y_final - 80 * (1.0 - anim_t))
        a       = int(255 * min(1.0, local_t * 4))
        _draw_read(img, REF_X, y, REF_W, READ_H, a, READ_MARKS[idx])

    draw = ImageDraw.Draw(img)
    la = int(255 * _ease_io(max(0.0, (t - 0.58) / 0.42)))
    _scene_text(draw, "Sequence",
                "Long-read whole-plasmid sequencing captures every gene copy", la)
    return img


# ═══════════════════════════════════════════════════════════════════════════════
# SCENE 2 — ANALYSE
# ═══════════════════════════════════════════════════════════════════════════════

def scene_analyse(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)

    _draw_ref_bar(img, 110, 38, 980, 18, 255)

    # Chart bounds
    CX, CY, CW, CH = 168, 90, 910, 325
    CB = CY + CH  # = 415

    roi_x1 = CX + int(CW * ROI_START)
    roi_x2 = CX + int(CW * ROI_END)

    ax_a = int(255 * min(1.0, t * 5))

    # ROI shading
    ov = Image.new("RGBA", img.size, (0,0,0,0))
    ImageDraw.Draw(ov).rectangle((roi_x1, CY, roi_x2, CB), fill=(*GREEN, 20))
    img.alpha_composite(ov)

    draw = ImageDraw.Draw(img)
    # Axes
    draw.line([(CX, CY-4), (CX, CB+2)], fill=(*MID, ax_a), width=2)
    draw.line([(CX-2, CB), (CX+CW, CB)], fill=(*MID, ax_a), width=2)
    # Tick marks and labels
    for rate, lbl in [(0.005, "0.005"), (0.010, "0.010"), (0.020, "0.020")]:
        ty = _profile_to_y(rate, CY, CB)
        draw.line([(CX-5, ty), (CX+3, ty)], fill=(*MID, ax_a), width=1)
        tw = _textw(F_SM, lbl)
        draw.text((CX - tw - 7, ty - 8), lbl, font=F_SM, fill=(*MID, ax_a))

    # Background rate band + label
    bg_y = _profile_to_y(BACKGROUND_RATE, CY, CB)
    bov  = Image.new("RGBA", img.size, (0,0,0,0))
    ImageDraw.Draw(bov).rectangle((CX, bg_y-8, CX+CW, bg_y+8), fill=(*SLATE, 28))
    img.alpha_composite(bov)
    draw = ImageDraw.Draw(img)
    draw.line([(CX, bg_y), (CX+CW, bg_y)], fill=(*SLATE, 140), width=1)
    draw.text((CX+CW+5, bg_y-8), "background", font=F_SM, fill=(*SLATE, ax_a))

    # Axis labels
    lbl = "mismatch rate"
    draw.text((CX - _textw(F_SM, lbl)//2 - 20, CY - 22), lbl, font=F_SM, fill=(*MID, ax_a))
    draw.text((roi_x1 + (roi_x2-roi_x1)//2 - 14, CY-18), "gene", font=F_SM, fill=(*GREEN, ax_a))
    draw.text((CX+CW//2 - 30, CB+10), "position (bp)", font=F_SM, fill=(*MID, ax_a))

    # Animated line trace
    line_t = _ease_out(max(0.0, (t - 0.04) / 0.96))
    n_pts  = int(PROFILE_LEN * line_t)
    if n_pts > 1:
        pts = []
        for i in range(n_pts):
            px = CX + int(CW * i / PROFILE_LEN)
            py = _profile_to_y(_PROFILE[i], CY, CB)
            pts.append((px, py))
        draw.line(pts, fill=(*BLUE, 210), width=2)

    la = int(255 * _ease_io(max(0.0, (t - 0.62) / 0.38)))
    _scene_text(draw, "Analyse",
                "Per-base mismatch rate reveals elevated mutation within the gene", la)
    return img


# ═══════════════════════════════════════════════════════════════════════════════
# SCENE 3 — MODEL
# ═══════════════════════════════════════════════════════════════════════════════

def scene_model(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)
    draw = ImageDraw.Draw(img)

    # Lambda counter
    la_main = int(255 * min(1.0, t * 3.5))
    lam_val  = LAMBDA_AA * _ease_out(min(1.0, t * 2.0))
    lam_str  = f"λ = {lam_val:.1f}"
    lam_x    = 90
    draw.text((lam_x, 50),  "rate  ×  gene length  →",   font=F_MED, fill=(*MID, la_main))
    draw.text((lam_x, 82),  lam_str,                      font=F_BIG, fill=(*INK, la_main))
    draw.text((lam_x, 138), "mutations per gene copy",    font=F_MED, fill=(*MID, la_main))

    # Separator
    draw.line([(560, 45), (560, 440)], fill=(*BORD, la_main), width=1)

    # Six example gene copies (3 per column, 2 rows) — right panel x=580..1110
    COPY_W, COPY_H = 238, 26
    COL_X = [585, 585 + COPY_W + 28]
    ROW_Y = [85, 175, 270, 85+90*3]

    layouts = [
        (0, 0), (1, 0),
        (0, 1), (1, 1),
        (0, 2), (1, 2),
    ]

    for ci, (col, row) in enumerate(layouts):
        cx = COL_X[col]
        cy = 70 + row * 110
        delay   = 0.04 * ci + 0.15
        local_t = _ease_out(max(0.0, (t - delay) / max(1e-4, 1.0 - delay)))
        a = int(255 * local_t)
        if a < 5:
            continue

        # Gene bar (green)
        draw.rounded_rectangle((cx, cy, cx+COPY_W, cy+COPY_H), radius=COPY_H//2,
                                fill=(*GREEN, a), outline=(*BORD, a//2), width=1)
        # Mutation dots
        for frac in _COPY_MUTATIONS[ci]:
            mx = cx + int(COPY_W * frac)
            my = cy + COPY_H // 2
            r  = 6
            draw.ellipse((mx-r, my-r, mx+r, my+r), fill=(*CORAL, a))
            draw.ellipse((mx-r+2, my-r+2, mx+r-2, my+r-2),
                         outline=(255,255,255,a), width=1)

        # Label beneath
        if local_t > 0.45:
            ta  = int(a * (local_t - 0.45) / 0.55)
            lbl = _COPY_LABELS[ci]
            tw  = _textw(F_SM, lbl)
            draw.text((cx + COPY_W//2 - tw//2, cy + COPY_H + 4),
                      lbl, font=F_SM, fill=(*MID, ta))

    la = int(255 * _ease_io(max(0.0, (t - 0.60) / 0.40)))
    _scene_text(draw, "Model",
                "Monte Carlo translation predicts the amino-acid mutation burden", la)
    return img


# ═══════════════════════════════════════════════════════════════════════════════
# SCENE 4 — REPORT
# ═══════════════════════════════════════════════════════════════════════════════

def scene_report(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)
    draw = ImageDraw.Draw(img)

    CX, CB, CW, CH = 185, 425, 840, 315
    CY     = CB - CH   # = 110
    MAX_V  = max(HIST_VALS) * 1.15
    BAR_W  = CW // 5 - 24   # ≈ 144 px
    GAP    = CW // 5         # ≈ 168 px per slot

    ax_a = int(255 * min(1.0, t * 5))
    draw.line([(CX, CY-4), (CX, CB+2)], fill=(*MID, ax_a), width=2)
    draw.line([(CX-2, CB), (CX+CW, CB)], fill=(*MID, ax_a), width=2)

    # Y-axis ticks
    for frac, lbl in [(0.1,"10 %"), (0.2,"20 %"), (0.3,"30 %")]:
        ty = int(CB - frac / MAX_V * CH)
        draw.line([(CX-5, ty), (CX+3, ty)], fill=(*MID, ax_a), width=1)
        tw = _textw(F_SM, lbl)
        draw.text((CX - tw - 7, ty - 8), lbl, font=F_SM, fill=(*MID, ax_a))

    draw.text((CX-10, CY-22), "fraction of library", font=F_SM, fill=(*MID, ax_a))
    xlab = "AA mutations per gene copy"
    draw.text((CX + CW//2 - _textw(F_SM,xlab)//2, CB+10), xlab,
              font=F_SM, fill=(*MID, ax_a))

    # Bars
    for k, (val, lbl, col) in enumerate(zip(HIST_VALS, HIST_LABELS, HIST_COLS)):
        delay = 0.03 * k
        bar_t = _ease_out(max(0.0, (t - delay - 0.05) / 0.95))
        bx    = CX + k * GAP + 12
        bar_h = int(CH * val / MAX_V * bar_t)

        if bar_h > 0:
            draw.rounded_rectangle((bx, CB-bar_h, bx+BAR_W, CB),
                                   radius=6, fill=(*col, 210))

        if bar_t > 0.65:
            ta      = int(210 * (bar_t - 0.65) / 0.35)
            pct_str = f"{val*100:.0f}%"
            tw      = _textw(F_SM, pct_str)
            draw.text((bx + BAR_W//2 - tw//2, CB - bar_h - 22),
                      pct_str, font=F_SM, fill=(*INK, ta))

        draw.text((bx + BAR_W//2 - _textw(F_SM,lbl)//2, CB+8),
                  lbl, font=F_SM, fill=(*MID, ax_a))

    # Lambda annotation — placed below the tallest bars
    if t > 0.68:
        ta  = int(255 * (t - 0.68) / 0.32)
        ann = f"λ = {LAMBDA_AA:.1f}   ·   {HIST_VALS[0]*100:.0f}% wild-type"
        draw.text((CX + CW//2 - _textw(F_MED, ann)//2, CB + 36), ann,
                  font=F_MED, fill=(*BLUE, ta))

    la = int(255 * _ease_io(max(0.0, (t - 0.58) / 0.42)))
    _scene_text(draw, "Report",
                "Poisson model reveals the fraction of wild-type and mutant variants", la)
    return img


# ═══════════════════════════════════════════════════════════════════════════════
# ASSEMBLE
# ═══════════════════════════════════════════════════════════════════════════════

def _crossfade(a: Image.Image, b: Image.Image, t: float) -> Image.Image:
    return Image.blend(a.convert("RGBA"), b.convert("RGBA"), _ease_io(t))


def build_frames() -> List[Image.Image]:
    fns      = [scene_sequence, scene_analyse, scene_model, scene_report]
    rendered = [[fn(i) for i in range(SCENE_FRAMES)] for fn in fns]
    frames: List[Image.Image] = []
    for si, sf in enumerate(rendered):
        frames.extend(sf)
        if si < len(rendered) - 1:
            frames.extend([sf[-1]] * HOLD_FRAMES)
            for fi in range(FADE_FRAMES):
                t = (fi + 1) / (FADE_FRAMES + 1)
                frames.append(_crossfade(sf[-1], rendered[si+1][0], t))
    return [f.convert("P", palette=Image.Palette.ADAPTIVE, colors=255)
            for f in frames]


_STATIC_COPY = Path("src/uht_tooling/web/static/animations/ep_library_profiler.gif")


def main() -> None:
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    frames = build_frames()
    scene_report(SCENE_FRAMES - 1).convert("RGB").save(POSTER_PATH)
    frames[0].save(
        OUTPUT_PATH,
        save_all=True,
        append_images=frames[1:],
        duration=int(1000 / FPS),
        loop=0,
        optimize=False,
        disposal=2,
    )
    print(f"✓  {len(frames)} frames → {OUTPUT_PATH}")
    if _STATIC_COPY.parent.exists():
        import shutil
        shutil.copy2(OUTPUT_PATH, _STATIC_COPY)
        print(f"✓  synced → {_STATIC_COPY}")


if __name__ == "__main__":
    main()
