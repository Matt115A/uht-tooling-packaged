"""Minimal, beautiful explainer GIF — Synthetic Gene Pool workflow.

Light-mode Apple style + Super Mario-ish bounce physics:
  • Bars drop 80 px and elastic-bounce on landing (spring overshoot)
  • Coin-collect sparkle on the unique-index barcode when each bar lands
  • Primer arrows shoot in and snap with a contact flash
  • Expanding spotlight ring on the selected bar
  • Protein dots pop out like coins collected
"""

from __future__ import annotations
import math
from pathlib import Path
from typing import List, Tuple

from PIL import Image, ImageDraw, ImageFilter, ImageFont

WIDTH        = 1200
HEIGHT       = 675
FPS          = 24
SCENE_FRAMES = 48
FADE_FRAMES  = 12
HOLD_FRAMES  = 6   # still frames held between scenes

OUTPUT_PATH = Path("assets/animations/synthetic_gene_pool_workflow.gif")
POSTER_PATH = Path("assets/animations/synthetic_gene_pool_workflow_poster.png")

# ── GUI colour tokens (light mode, mirrors theme.py) ─────────────────────────
BG     = (245, 245, 247)   # --bg-secondary  #f5f5f7
CARD   = (255, 255, 255)   # --bg-primary    #ffffff
INK    = (29,  29,  31)    # --text-primary  #1d1d1f
MID    = (110, 110, 115)   # --text-secondary #6e6e73
BLUE   = (0,   113, 227)   # --accent        #0071e3
BORD   = (210, 210, 215)   # card border
BAR_BG = (248, 248, 250)   # inner-bar background

# ── DNA segment colours ───────────────────────────────────────────────────────
TEAL   = (0,   196, 175)
GREEN  = (40,  185, 115)
AMBER  = (255, 180, 50)
CORAL  = (255, 80,  65)
VIOLET = (155, 100, 240)
SLATE  = (75,  115, 160)

# ── Construct layout ──────────────────────────────────────────────────────────
# Segments: [5′-handle, gene, RBS, spacer, barcode, 3′-handle]
SEG_FRACS  = [0.18, 0.28, 0.10, 0.07, 0.15, 0.22]
SEG_COLORS = [TEAL, GREEN, AMBER, CORAL, VIOLET, SLATE]
NUM_TARGETS = 6
SELECT_IDX  = 2

# Cumulative fraction to centre of barcode (VIOLET) segment
BARCODE_CUM = 0.18 + 0.28 + 0.10 + 0.07 + 0.075   # = 0.705

# ── Per-member barcode patterns ───────────────────────────────────────────────
# Each row = alternating (dark, light) stripe widths in relative units.
# First stripe is always dark. Total = 24 units each — gives distinct patterns
# without looking cluttered.
_BARCODE_PATS: List[List[int]] = [
    [3, 2, 1, 1, 4, 2, 1, 2, 3, 2, 3],   # member 0
    [1, 1, 4, 1, 2, 3, 2, 1, 4, 2, 3],   # member 1
    [2, 2, 2, 1, 1, 2, 4, 1, 3, 2, 4],   # member 2
    [4, 1, 2, 2, 1, 3, 2, 1, 2, 2, 4],   # member 3
    [1, 3, 3, 1, 2, 2, 1, 1, 4, 2, 4],   # member 4
    [3, 1, 1, 2, 3, 1, 1, 3, 3, 2, 4],   # member 5
]
_BARCODE_TOTAL = 24   # sum of each row above

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
    """Spring-bounce: overshoots ~18 % then settles. Mario-landing feel."""
    if t <= 0: return 0.0
    if t >= 1: return 1.0
    p, s = 0.38, 0.095
    return 2 ** (-9 * t) * math.sin((t - s) * 2 * math.pi / p) + 1


def _shoot_ease(t: float) -> float:
    """Fast cubic approach + tiny elastic snap at the end."""
    t    = max(0.0, min(1.0, t))
    base = 1 - (1 - t) ** 3
    snap = 0.09 * math.exp(-22 * (t - 0.82) ** 2) if t > 0.70 else 0.0
    return base + snap

# ── Primitive helpers ─────────────────────────────────────────────────────────

def _base_img() -> Image.Image:
    return Image.new("RGBA", (WIDTH, HEIGHT), BG)


def _draw_card(img: Image.Image, x0: int, y0: int,
               x1: int, y1: int, radius: int = 18) -> None:
    sh = Image.new("RGBA", img.size, (0, 0, 0, 0))
    sd = ImageDraw.Draw(sh)
    sd.rounded_rectangle((x0 + 4, y0 + 8, x1 + 4, y1 + 8),
                          radius=radius, fill=(0, 0, 0, 32))
    img.alpha_composite(sh.filter(ImageFilter.GaussianBlur(12)))
    ImageDraw.Draw(img).rounded_rectangle(
        (x0, y0, x1, y1), radius=radius,
        fill=(*CARD, 255), outline=(*BORD, 255), width=1)


def _draw_barcode_stripes(draw: ImageDraw.ImageDraw,
                           zx: int, zy: int, zw: int, zh: int,
                           member_idx: int, alpha: int) -> None:
    """Draw vertical barcode stripes for a specific pool member."""
    pat  = _BARCODE_PATS[member_idx % len(_BARCODE_PATS)]
    unit = zw / _BARCODE_TOTAL
    margin = 5          # vertical inset so stripes don't bleed past rounded ends
    cursor = zx
    is_dark = True      # first stripe is always dark
    for width_units in pat:
        w = unit * width_units
        if is_dark:
            x0 = int(cursor) + 1
            x1 = int(cursor + w) - 1
            draw.rectangle((x0, zy + margin, max(x0 + 1, x1), zy + zh - margin),
                            fill=(*VIOLET, alpha))
        cursor += w
        is_dark = not is_dark


def _draw_construct(img: Image.Image, x: int, y: int, w: int, h: int,
                    alpha: int = 255, highlight: bool = False,
                    progress: float = 1.0, member_idx: int = 0) -> None:
    pw = max(0, min(w, int(w * progress)))
    if pw < 6:
        return
    draw = ImageDraw.Draw(img)
    draw.rounded_rectangle((x, y, x + pw, y + h), radius=h // 2,
                            fill=(*BAR_BG, alpha),
                            outline=(*BORD, alpha), width=1)
    cursor, inner = x + 4, pw - 8
    gene_bounds: Tuple[int, int] | None = None

    for i, (frac, col) in enumerate(zip(SEG_FRACS, SEG_COLORS)):
        sw = int(inner * frac)
        if i == len(SEG_FRACS) - 1:
            sw = (x + pw - 4) - cursor
        if sw < 2:
            cursor += sw + 2
            continue

        if i == 4:  # barcode segment — unique per member
            # Light background with small radius so stripes don't bleed
            draw.rounded_rectangle((cursor, y + 4, cursor + sw, y + h - 4),
                                    radius=5, fill=(*BAR_BG, alpha))
            _draw_barcode_stripes(draw, cursor, y + 4, sw, h - 8,
                                   member_idx, alpha)
        else:
            draw.rounded_rectangle((cursor, y + 4, cursor + sw, y + h - 4),
                                    radius=(h - 8) // 2, fill=(*col, alpha))

        if i == 1:  # GREEN = gene segment
            gene_bounds = (cursor, cursor + sw)

        cursor += sw + 2

    # "gene" label centred inside the gene segment
    if gene_bounds and alpha > 60:
        gx1, gx2 = gene_bounds
        gcx = (gx1 + gx2) // 2
        label = "gene"
        try:
            tw = int(F_TINY.getlength(label))
        except AttributeError:
            tw = 28
        ty_label = y + h // 2 - 6
        label_alpha = min(255, int(alpha * 0.95))
        draw.text((gcx - tw // 2, ty_label), label,
                  font=F_TINY, fill=(255, 255, 255, label_alpha))

    if highlight:
        g  = Image.new("RGBA", img.size, (0, 0, 0, 0))
        gd = ImageDraw.Draw(g)
        gd.rounded_rectangle((x - 10, y - 10, x + pw + 10, y + h + 10),
                              radius=h // 2 + 10, fill=(*BLUE, 20))
        img.alpha_composite(g.filter(ImageFilter.GaussianBlur(14)))
        ImageDraw.Draw(img).rounded_rectangle(
            (x, y, x + pw, y + h), radius=h // 2,
            outline=(*BLUE, min(255, alpha)), width=3)


def _arrow(draw: ImageDraw.ImageDraw, x0: float, y0: float,
           x1: float, y1: float, color: Tuple, width: int = 5) -> None:
    draw.line(((x0, y0), (x1, y1)), fill=color, width=width)
    ang  = math.atan2(y1 - y0, x1 - x0)
    sz   = 13
    tip  = (x1, y1)
    draw.polygon([tip,
                  (tip[0] - sz * math.cos(ang - math.pi / 7),
                   tip[1] - sz * math.sin(ang - math.pi / 7)),
                  (tip[0] - sz * math.cos(ang + math.pi / 7),
                   tip[1] - sz * math.sin(ang + math.pi / 7))],
                 fill=color)


def _sparkles(img: Image.Image, cx: int, cy: int, t: float,
              color: Tuple, n: int = 8, dist: int = 30) -> None:
    draw = ImageDraw.Draw(img)
    for i in range(n):
        ang = i * 2 * math.pi / n + math.pi / n
        d   = _ease_out(t) * dist
        r   = max(1, int(5 * (1 - t)))
        a   = int(230 * (1 - t))
        px, py = int(cx + d * math.cos(ang)), int(cy + d * math.sin(ang))
        draw.ellipse((px - r, py - r, px + r, py + r), fill=(*color, a))


def _scene_text(draw: ImageDraw.ImageDraw,
                title: str, subtitle: str, alpha: int) -> None:
    draw.rectangle((70, 488, WIDTH - 70, 490), fill=(*BLUE, alpha // 3))
    draw.text((80, 502),  title,    font=F_HUGE, fill=(*INK, alpha))
    draw.text((82, 588),  subtitle, font=F_MED,  fill=(*MID, alpha))


def _pool() -> Tuple[int, int, int, int, int]:
    """sx, ty, bw, bh, spacing"""
    return 110, 76, 980, 44, 62


# ═══════════════════════════════════════════════════════════════════════════════
# SCENE 1 — DESIGN
# ═══════════════════════════════════════════════════════════════════════════════

def scene_design(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)
    sx, ty, bw, bh, sp = _pool()

    for idx in range(NUM_TARGETS):
        delay   = idx * 0.072
        local_t = max(0.0, min(1.0, (t - delay) / max(1e-4, 1.0 - delay)))
        anim_t  = _elastic_out(local_t)

        y_final = ty + idx * sp
        y = int(y_final - 80 * (1.0 - anim_t))
        a = int(255 * min(1.0, local_t * 3.5))
        _draw_construct(img, sx, y, bw, bh, alpha=a, member_idx=idx)

        # Sparkle at barcode segment when bar lands
        if local_t > 0.58:
            spark_t = _ease_out(min(1.0, (local_t - 0.58) / 0.42))
            vx = sx + 4 + int((bw - 8) * BARCODE_CUM)
            _sparkles(img, vx, y_final + bh // 2, spark_t, VIOLET)

    draw = ImageDraw.Draw(img)
    la = int(255 * _ease_io(max(0.0, (t - 0.54) / 0.46)))
    _scene_text(draw, "Design",
                "Pool of gene constructs — each with a unique barcode index", la)
    return img


# ═══════════════════════════════════════════════════════════════════════════════
# SCENE 2 — AMPLIFY
# ═══════════════════════════════════════════════════════════════════════════════

def scene_amplify(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)
    sx, ty, bw, bh, sp = _pool()

    for idx in range(NUM_TARGETS):
        _draw_construct(img, sx, ty + idx * sp, bw, bh, member_idx=idx)

    pt = _ease_out(max(0.0, (t - 0.18) / 0.82))

    fw_x1  = int(sx + bw * 0.18 * _shoot_ease(pt))
    rv_x0  = int(sx + bw - bw * 0.22 * _shoot_ease(pt))

    draw = ImageDraw.Draw(img)
    if pt > 0.02:
        fy = ty - 38
        _arrow(draw, sx - 10, fy, fw_x1, fy, (*BLUE, int(240 * pt)))
        if pt > 0.68:
            ta = int(210 * (pt - 0.68) / 0.32)
            draw.text((sx, fy - 22), "POOL_CONST_F  →", font=F_SM,
                      fill=(*BLUE, ta))

        ry = min(ty + NUM_TARGETS * sp + 16, 455)
        _arrow(draw, sx + bw + 10, ry, rv_x0, ry, (*SLATE, int(240 * pt)))
        if pt > 0.68:
            ta = int(210 * (pt - 0.68) / 0.32)
            draw.text((rv_x0 - 4, ry + 9), "←  POOL_CONST_R",
                      font=F_SM, fill=(*SLATE, ta))

        if 0.84 < pt < 0.99:
            flash_t = (pt - 0.84) / 0.15
            fa = int(150 * (1 - flash_t))
            fr = max(2, int(20 * (1 - flash_t)))
            fl = Image.new("RGBA", img.size, (0, 0, 0, 0))
            fd = ImageDraw.Draw(fl)
            fd.ellipse((fw_x1 - fr, fy - fr, fw_x1 + fr, fy + fr),
                        fill=(*BLUE, fa))
            fd.ellipse((rv_x0 - fr, ry - fr, rv_x0 + fr, ry + fr),
                        fill=(*SLATE, fa))
            img.alpha_composite(fl.filter(ImageFilter.GaussianBlur(5)))
            draw = ImageDraw.Draw(img)

    la = int(255 * _ease_io(max(0.0, (t - 0.54) / 0.46)))
    _scene_text(draw, "Amplify",
                "One shared primer pair amplifies the entire pool", la)
    return img


# ═══════════════════════════════════════════════════════════════════════════════
# SCENE 3 — SELECT
# ═══════════════════════════════════════════════════════════════════════════════

def scene_select(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)
    sx, ty, bw, bh, sp = _pool()

    for idx in range(NUM_TARGETS):
        _draw_construct(img, sx, ty + idx * sp, bw, bh, member_idx=idx)

    dim_t = _ease_io(min(1.0, t * 1.8))
    if dim_t > 0.01:
        ov = Image.new("RGBA", img.size, (0, 0, 0, 0))
        od = ImageDraw.Draw(ov)
        for idx in range(NUM_TARGETS):
            if idx != SELECT_IDX:
                y = ty + idx * sp
                od.rounded_rectangle((sx, y, sx + bw, y + bh),
                                     radius=bh // 2,
                                     fill=(*BG, int(215 * dim_t)))
        img.alpha_composite(ov)

    _draw_construct(img, sx, ty + SELECT_IDX * sp, bw, bh,
                    highlight=t > 0.26, member_idx=SELECT_IDX)
    draw = ImageDraw.Draw(img)

    if 0.20 < t < 0.88:
        ring_t = _ease_out(min(1.0, (t - 0.20) / 0.46))
        cx = sx + bw // 2
        cy = ty + SELECT_IDX * sp + bh // 2
        r  = int(48 + ring_t * 105)
        ra = int(72 * (1 - ring_t))
        ri = Image.new("RGBA", img.size, (0, 0, 0, 0))
        ImageDraw.Draw(ri).ellipse((cx - r, cy - r, cx + r, cy + r),
                                   outline=(*BLUE, ra), width=3)
        img.alpha_composite(ri)
        draw = ImageDraw.Draw(img)

    # ── Pullout primer ────────────────────────────────────────────────────────
    # Arrow runs ABOVE the selected bar so it's never hidden by the construct.
    # A short vertical tick drops from the arrow tip to the bar's barcode zone.
    pull_raw = max(0.0, (t - 0.26) / 0.74)
    pt = _elastic_out(pull_raw)
    if pull_raw > 0.02:
        bar_top  = ty + SELECT_IDX * sp          # top edge of selected bar
        primer_y = bar_top - 20                  # arrow lane: 20 px above bar

        # Barcode left-edge x (VIOLET/barcode starts at ~70.5 % from bar left)
        barcode_x = sx + 4 + int((bw - 8) * (BARCODE_CUM - 0.075))

        pull_x0 = sx + bw + 10 + int(88 * max(0.0, 1 - pull_raw * 1.5))
        pull_x1 = int(sx + bw - bw * 0.30 * pt)

        primer_alpha = int(240 * min(1.0, pull_raw * 1.6))

        if pull_x1 < pull_x0 - 5:
            # Horizontal approach above the bar
            _arrow(draw, pull_x0, primer_y, pull_x1, primer_y,
                   (*VIOLET, primer_alpha))

        # Downward tick from arrow tip to bar top — appears once primer is close
        if pull_raw > 0.72 and pull_x1 < pull_x0 - 5:
            tick_a = int(primer_alpha * min(1.0, (pull_raw - 0.72) / 0.18))
            draw.line([(pull_x1, primer_y), (pull_x1, bar_top)],
                      fill=(*VIOLET, tick_a), width=3)
            # Small dot where tick meets the bar top
            r = 4
            draw.ellipse((pull_x1 - r, bar_top - r, pull_x1 + r, bar_top + r),
                          fill=(*VIOLET, tick_a))

        if t > 0.68:
            ta = int(220 * (t - 0.68) / 0.32)
            draw.text((sx + bw + 14, primer_y - 20), "pullout primer",
                      font=F_SM, fill=(*VIOLET, ta))

    la = int(255 * _ease_io(max(0.0, (t - 0.54) / 0.46)))
    _scene_text(draw, "Select",
                "A gene-specific primer recovers one variant from the pool", la)
    return img


# ═══════════════════════════════════════════════════════════════════════════════
# SCENE 4 — EXPRESS
# ═══════════════════════════════════════════════════════════════════════════════

def scene_express(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)

    bw, bh = 820, 50
    bx     = (WIDTH - bw) // 2
    by     = 96
    _draw_construct(img, bx, by, bw, bh, highlight=True, member_idx=SELECT_IDX)

    draw = ImageDraw.Draw(img)
    draw.text((bx, by + bh + 10), "CFPS-ready — no cloning required",
              font=F_SM, fill=(*TEAL, 210))

    ax  = bx + bw // 2
    ay0 = by + bh + 6
    ay1 = int(ay0 + 74 * _ease_out(t))
    if ay1 > ay0 + 5:
        _arrow(draw, ax, ay0, ax, ay1, (*GREEN, 205))

    ft = _ease_out(max(0.0, (t - 0.12) / 0.88))
    if ft > 0.04:
        fx, fw = ax - 40, 80
        fy     = ay0 + 80
        fh     = int(144 * ft)

        sh = Image.new("RGBA", img.size, (0, 0, 0, 0))
        sd = ImageDraw.Draw(sh)
        sd.rounded_rectangle((fx + 4, fy + 8, fx + fw + 4, fy + fh + 8),
                              radius=16, fill=(0, 0, 0, 24))
        img.alpha_composite(sh.filter(ImageFilter.GaussianBlur(9)))

        draw = ImageDraw.Draw(img)
        draw.rounded_rectangle((fx, fy, fx + fw, fy + fh), radius=16,
                                fill=(*CARD, 252), outline=(*BORD, 255), width=2)
        fill_h = int(fh * 0.50 * ft)
        if fill_h > 4:
            draw.rectangle((fx + 5, fy + fh - fill_h,
                             fx + fw - 5, fy + fh - 5),
                            fill=(*GREEN, 130))
        draw.polygon([(fx + 14, fy + fh), (fx + fw - 14, fy + fh),
                      (fx + fw // 2, fy + fh + 28)],
                     fill=(*CARD, 252), outline=(*BORD, 255))

    COIN_COLS = [BLUE, GREEN, VIOLET, AMBER, TEAL, CORAL]
    center_y  = ay0 + 158
    for i in range(12):
        prog = _ease_out(max(0.0, (t - 0.36 - i * 0.018) / 0.64))
        if prog > 0:
            ang = i * 2 * math.pi / 12
            px  = int(ax + 95 * math.cos(ang) * prog)
            py  = int(center_y + 68 * math.sin(ang) * prog)
            r   = max(1, int((7 + i % 3) * (1 - prog * 0.45)))
            a   = int(235 * (1 - prog * 0.40))
            draw.ellipse((px - r, py - r, px + r, py + r),
                          fill=(*COIN_COLS[i % len(COIN_COLS)], a))

    la = int(255 * _ease_io(max(0.0, (t - 0.54) / 0.46)))
    _scene_text(draw, "Express",
                "Recovered product goes straight into the cell-free reaction", la)
    return img


# ═══════════════════════════════════════════════════════════════════════════════
# ASSEMBLE
# ═══════════════════════════════════════════════════════════════════════════════

def _crossfade(a: Image.Image, b: Image.Image, t: float) -> Image.Image:
    return Image.blend(a.convert("RGBA"), b.convert("RGBA"), _ease_io(t))


def build_frames() -> List[Image.Image]:
    fns      = [scene_design, scene_amplify, scene_select, scene_express]
    rendered = [[fn(i) for i in range(SCENE_FRAMES)] for fn in fns]
    frames: List[Image.Image] = []
    for si, sf in enumerate(rendered):
        frames.extend(sf)
        if si < len(rendered) - 1:
            # Hold on last frame before crossfading
            frames.extend([sf[-1]] * HOLD_FRAMES)
            for fi in range(FADE_FRAMES):
                t = (fi + 1) / (FADE_FRAMES + 1)
                frames.append(_crossfade(sf[-1], rendered[si + 1][0], t))
    return [f.convert("P", palette=Image.Palette.ADAPTIVE, colors=255)
            for f in frames]


_STATIC_COPY = Path("src/uht_tooling/web/static/animations/synthetic_gene_pool_workflow.gif")


def main() -> None:
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    frames = build_frames()
    scene_express(SCENE_FRAMES - 1).convert("RGB").save(POSTER_PATH)
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
