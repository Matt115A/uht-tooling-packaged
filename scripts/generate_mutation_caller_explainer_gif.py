"""Mutation Caller — explainer GIF.

Scenes:
  1. Locate   — reads drop in; upstream/downstream flanks light up; gene region is boxed
  2. Extract  — gene regions lift out from each read and stack below
  3. Align    — extracted copies align to the template; mutation marks appear
  4. Report   — AA substitution frequency bar chart grows up
"""

from __future__ import annotations
import math
import random
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from PIL import Image, ImageDraw, ImageFilter, ImageFont

WIDTH, HEIGHT = 1200, 675
FPS          = 24
SCENE_FRAMES = 48
FADE_FRAMES  = 12
HOLD_FRAMES  = 6

OUTPUT_PATH = Path("assets/animations/mutation_caller.gif")
POSTER_PATH = Path("assets/animations/mutation_caller_poster.png")

# ── Colours ───────────────────────────────────────────────────────────────────
BG     = (245, 245, 247)
CARD   = (255, 255, 255)
INK    = (29,  29,  31)
MID    = (110, 110, 115)
BLUE   = (0,   113, 227)
BORD   = (210, 210, 215)
TEAL   = (0,   196, 175)
GREEN  = (40,  185, 115)
AMBER  = (255, 180, 50)
CORAL  = (255, 80,  65)
VIOLET = (155, 100, 240)
SLATE  = (75,  115, 160)

# ── Read anatomy fractions (proportions of bar inner width) ───────────────────
# [vector left | upstream flank | gene | downstream flank | vector right]
_VL  = 0.14   # vector left width
_UF  = 0.08   # upstream flank width
_GN  = 0.52   # gene width
_DF  = 0.08   # downstream flank width
_VR  = 0.18   # vector right width   (total = 1.0)

UPS_START  = _VL
UPS_END    = _VL + _UF
GENE_START = UPS_END
GENE_END   = UPS_END + _GN
DWN_START  = GENE_END
DWN_END    = GENE_END + _DF

# ── Scene 3 mutation positions (fraction along gene bar) ──────────────────────
# Three recurring sites — maps to A42V, G78R, L15P in the report
MUT_FRAC = [0.12, 0.35, 0.63]

# Which reads carry which mutations (5 reads in align scene)
ALIGN_MUTS: List[List[int]] = [
    [0, 1],     # read 0: sites 0 and 1
    [1],        # read 1: site 1 only
    [],         # read 2: WT
    [0],        # read 3: site 0 only
    [2],        # read 4: site 2 only
]

# ── Scene 4 bar chart data ────────────────────────────────────────────────────
MC_LABELS = ["A42V", "G78R", "L15P", "D91G", "V53I", "R120K", "T8A"]
MC_COUNTS = [245,    189,    156,    98,     87,     42,      18]
MC_THRESH = 50   # minimum count to be "frequent" (highlighted)

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
F_BIG  = _font(36, bold=True)
F_MED  = _font(20)
F_SM   = _font(14)
F_TINY = _font(11, bold=True)
F_XS   = _font(10)

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
    sh = Image.new("RGBA", img.size, (0,0,0,0))
    sd = ImageDraw.Draw(sh)
    sd.rounded_rectangle((x0+4, y0+8, x1+4, y1+8), radius=18, fill=(0,0,0,32))
    img.alpha_composite(sh.filter(ImageFilter.GaussianBlur(12)))
    ImageDraw.Draw(img).rounded_rectangle(
        (x0, y0, x1, y1), radius=18,
        fill=(*CARD,255), outline=(*BORD,255), width=1)

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


def _seg_xs(x: int, w: int) -> Dict[str, Tuple[int, int]]:
    """Pixel (left, right) for each anatomy segment given bar x and w."""
    inner = w - 8
    segs: Dict[str, Tuple[float, float]] = {
        "vl":   (0.0,       UPS_START),
        "up":   (UPS_START, UPS_END),
        "gene": (GENE_START, GENE_END),
        "dn":   (DWN_START, DWN_END),
        "vr":   (DWN_END,   1.0),
    }
    return {
        name: (x + 4 + int(inner * s), x + 4 + int(inner * e))
        for name, (s, e) in segs.items()
    }


def _draw_read_bar(img: Image.Image, x: int, y: int, w: int, h: int,
                   alpha: int = 255,
                   flank_glow: float = 0.0,
                   gene_highlight: bool = False) -> None:
    """Draw a full sequencing read showing its anatomy."""
    draw = ImageDraw.Draw(img)
    draw.rounded_rectangle((x, y, x+w, y+h), radius=h//2,
                            fill=(215,222,232,alpha), outline=(*BORD,alpha), width=1)
    sx = _seg_xs(x, w)
    inner_h = h - 8
    r = inner_h // 2

    def seg(name: str, col: Tuple[int,int,int]) -> None:
        x0, x1 = sx[name]
        if x1 > x0:
            draw.rounded_rectangle((x0, y+4, x1, y+h-4), radius=r,
                                   fill=(*col, alpha))

    seg("vl", SLATE)
    seg("dn", AMBER)
    seg("vr", SLATE)
    seg("gene", GREEN)

    # Gene label
    gx0, gx1 = sx["gene"]
    lbl = "gene"
    tw  = _textw(F_TINY, lbl)
    draw.text((gx0 + (gx1-gx0)//2 - tw//2, y+h//2-6), lbl,
              font=F_TINY, fill=(255,255,255,alpha))

    # Upstream flank drawn last so it appears on top
    ux0, ux1 = sx["up"]
    if ux1 > ux0:
        draw.rounded_rectangle((ux0, y+4, ux1, y+h-4), radius=r,
                               fill=(*AMBER, alpha))

    # Flank glow
    if flank_glow > 0.01:
        ga = int(80 * flank_glow)
        for seg_name in ("up", "dn"):
            fx0, fx1 = sx[seg_name]
            ov = Image.new("RGBA", img.size, (0,0,0,0))
            ImageDraw.Draw(ov).rounded_rectangle(
                (fx0-6, y-6, fx1+6, y+h+6), radius=r+6, fill=(*AMBER, ga))
            img.alpha_composite(ov.filter(ImageFilter.GaussianBlur(8)))

    # Gene highlight (blue border)
    if gene_highlight:
        gx0, gx1 = sx["gene"]
        ov = Image.new("RGBA", img.size, (0,0,0,0))
        ImageDraw.Draw(ov).rounded_rectangle(
            (gx0-3, y+1, gx1+3, y+h-1), radius=r+3, fill=(*BLUE, 18))
        img.alpha_composite(ov.filter(ImageFilter.GaussianBlur(6)))
        ImageDraw.Draw(img).rounded_rectangle(
            (gx0-3, y+1, gx1+3, y+h-1), radius=r+3,
            outline=(*BLUE, min(255,alpha)), width=2)


def _draw_gene_strip(img: Image.Image, x: int, y: int, w: int, h: int,
                     alpha: int = 255,
                     mut_sites: Optional[List[int]] = None) -> None:
    """Draw just the extracted gene (GREEN) with optional mutation marks."""
    draw = ImageDraw.Draw(img)
    draw.rounded_rectangle((x, y, x+w, y+h), radius=h//2,
                            fill=(*GREEN, alpha), outline=(*BORD, alpha//2), width=1)
    if mut_sites:
        for site_i in mut_sites:
            mx = x + int(w * MUT_FRAC[site_i])
            my = y + h // 2
            r  = 6
            draw.ellipse((mx-r, my-r, mx+r, my+r), fill=(*CORAL, alpha))
            draw.ellipse((mx-r+2, my-r+2, mx+r-2, my+r-2),
                         outline=(255,255,255,alpha), width=1)


# ═══════════════════════════════════════════════════════════════════════════════
# SCENE 1 — LOCATE
# ═══════════════════════════════════════════════════════════════════════════════

def scene_locate(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)

    NUM  = 6
    RX, RW, RH, RSP = 110, 980, 28, 54

    flank_glow = _ease_out(max(0.0, (t - 0.55) / 0.45))
    gene_hl    = t > 0.70

    for idx in range(NUM):
        delay   = idx * 0.08
        local_t = max(0.0, min(1.0, (t - delay) / max(1e-4, 1.0 - delay)))
        anim_t  = _elastic_out(local_t)
        y_final = 58 + idx * RSP
        y       = int(y_final - 80 * (1.0 - anim_t))
        a       = int(255 * min(1.0, local_t * 4))
        _draw_read_bar(img, RX, y, RW, RH, alpha=a,
                       flank_glow=flank_glow * min(1.0, local_t * 2),
                       gene_highlight=gene_hl)

    draw = ImageDraw.Draw(img)

    # Flank labels above first read (once reads have landed)
    if t > 0.60:
        la  = int(220 * (t - 0.60) / 0.40)
        sx  = _seg_xs(RX, RW)
        y0  = 58   # first read y

        # Upstream flank label
        ux0, ux1 = sx["up"]
        ucx = (ux0 + ux1) // 2
        draw.line([(ucx, y0 - 4), (ucx, y0 - 20)], fill=(*AMBER, la), width=2)
        lbl = "upstream flank"
        tw  = _textw(F_SM, lbl)
        draw.text((ucx - tw//2, y0 - 36), lbl, font=F_SM, fill=(*AMBER, la))

        # Downstream flank label
        dx0, dx1 = sx["dn"]
        dcx = (dx0 + dx1) // 2
        draw.line([(dcx, y0 - 4), (dcx, y0 - 20)], fill=(*AMBER, la), width=2)
        lbl2 = "downstream flank"
        tw2  = _textw(F_SM, lbl2)
        draw.text((dcx - tw2//2, y0 - 36), lbl2, font=F_SM, fill=(*AMBER, la))

        # Gene bracket
        gx0, gx1 = sx["gene"]
        if gene_hl:
            draw.line([(gx0, y0 - 4), (gx0, y0 - 16)], fill=(*BLUE, la), width=2)
            draw.line([(gx1, y0 - 4), (gx1, y0 - 16)], fill=(*BLUE, la), width=2)
            draw.line([(gx0, y0 - 16), (gx1, y0 - 16)], fill=(*BLUE, la), width=2)
            gcx = (gx0 + gx1) // 2
            draw.text((gcx - _textw(F_SM,"extracted region")//2, y0-32),
                      "extracted region", font=F_SM, fill=(*BLUE, la))

    la = int(255 * _ease_io(max(0.0, (t - 0.58) / 0.42)))
    _scene_text(draw, "Locate",
                "Flanking sequences anchor the gene inside every read", la)
    return img


# ═══════════════════════════════════════════════════════════════════════════════
# SCENE 2 — EXTRACT
# ═══════════════════════════════════════════════════════════════════════════════

def scene_extract(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)

    NUM   = 6
    RX, RW, RH, RSP = 110, 980, 22, 44

    # Top: static reads (dimmed)
    for idx in range(NUM):
        y = 42 + idx * RSP
        a = int(130 * (1.0 - min(1.0, t * 1.5)))
        if a > 5:
            _draw_read_bar(img, RX, y, RW, RH, alpha=a)

    # Gene strips lift out and fall into stack
    STACK_X  = 200
    STACK_W  = RW - 290
    STACK_H  = 26
    STACK_SP = 32

    draw = ImageDraw.Draw(img)

    for idx in range(NUM):
        delay   = idx * 0.07
        local_t = _ease_out(max(0.0, (t - delay - 0.05) / max(1e-4, 0.95 - delay)))
        a       = int(255 * min(1.0, local_t * 3))
        if a < 5:
            continue

        # Source y (from read position above)
        src_y = 42 + idx * RSP + (RH - STACK_H) // 2
        # Destination y (stacked in lower half of card)
        dst_y = 240 + idx * STACK_SP
        # Interpolated position
        y_cur = int(src_y + (dst_y - src_y) * local_t)

        # Source x (gene segment within full read)
        sx_map = _seg_xs(RX, RW)
        gx0, gx1 = sx_map["gene"]
        src_x = gx0
        # Destination x (centred stack)
        dst_x = STACK_X

        x_cur = int(src_x + (dst_x - src_x) * local_t)
        w_cur = STACK_W if local_t > 0.8 else int(gx1 - gx0)

        _draw_gene_strip(img, x_cur, y_cur, w_cur, STACK_H, alpha=a)

    # Label for the stack
    if t > 0.55:
        la = int(220 * (t - 0.55) / 0.45)
        draw = ImageDraw.Draw(img)
        lbl = "extracted gene copies"
        tw  = _textw(F_MED, lbl)
        draw.text((STACK_X + STACK_W + 16, 240 + NUM * STACK_SP // 2 - 10),
                  lbl, font=F_MED, fill=(*GREEN, la))

    draw = ImageDraw.Draw(img)
    la = int(255 * _ease_io(max(0.0, (t - 0.58) / 0.42)))
    _scene_text(draw, "Extract",
                "Gene copies are pulled from every read that spans both flanks", la)
    return img


# ═══════════════════════════════════════════════════════════════════════════════
# SCENE 3 — ALIGN
# ═══════════════════════════════════════════════════════════════════════════════

def scene_align(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)

    NUM_READS = len(ALIGN_MUTS)
    AX, AW, AH, ASP = 200, 800, 28, 56

    draw = ImageDraw.Draw(img)

    # Template (BLUE) at top
    tmpl_a = int(255 * min(1.0, t * 6))
    draw.rounded_rectangle((AX, 52, AX+AW, 52+AH), radius=AH//2,
                            fill=(*BLUE, tmpl_a), outline=(*BORD, tmpl_a//2), width=1)
    if tmpl_a > 80:
        lbl = "template"
        tw  = _textw(F_SM, lbl)
        draw.text((AX - tw - 10, 52 + AH//2 - 8), lbl,
                  font=F_SM, fill=(*BLUE, tmpl_a))

    # Reads slide in from the right and align
    for ri, muts in enumerate(ALIGN_MUTS):
        delay   = ri * 0.07 + 0.08
        local_t = _ease_out(max(0.0, (t - delay) / max(1e-4, 1.0 - delay)))
        a       = int(255 * min(1.0, local_t * 4))
        if a < 5:
            continue

        y_final = 52 + AH + 20 + ri * ASP
        x_final = AX
        x_start = AX + 320
        x_cur   = int(x_start + (x_final - x_start) * local_t)

        # Draw as a green gene strip (no mutations yet)
        _draw_gene_strip(img, x_cur, y_final, AW, AH, alpha=a, mut_sites=None)

        # Read label
        if local_t > 0.5:
            la2 = int(a * (local_t - 0.5) / 0.5)
            draw = ImageDraw.Draw(img)
            draw.text((x_cur - _textw(F_SM, f"read {ri+1}") - 10, y_final + AH//2 - 8),
                      f"read {ri+1}", font=F_SM, fill=(*MID, la2))

        # Mutation marks appear after read has settled
        if local_t > 0.82 and muts:
            mut_t = min(1.0, (local_t - 0.82) / 0.18)
            ma    = int(a * mut_t)
            draw = ImageDraw.Draw(img)
            for site_i in muts:
                mx = AX + int(AW * MUT_FRAC[site_i])
                my = y_final + AH // 2
                r  = int(6 * mut_t)
                if r > 0:
                    draw.ellipse((mx-r, my-r, mx+r, my+r), fill=(*CORAL, ma))
                if r >= 3:
                    draw.ellipse((mx-r+2, my-r+2, mx+r-2, my+r-2),
                                 outline=(255,255,255,ma), width=1)

    # Mutation site callouts (appear near end)
    if t > 0.78:
        ca = int(200 * (t - 0.78) / 0.22)
        draw = ImageDraw.Draw(img)
        for site_i, frac in enumerate(MUT_FRAC):
            mx = AX + int(AW * frac)
            # Vertical dashed indicator above template
            draw.line([(mx, 50), (mx, 42)], fill=(*CORAL, ca), width=2)
            # AA label
            labels_aa = ["L15P", "A42V", "G78R"]
            lbl = labels_aa[site_i]
            tw  = _textw(F_XS, lbl)
            draw.text((mx - tw//2, 30), lbl, font=F_XS, fill=(*CORAL, ca))

    draw = ImageDraw.Draw(img)
    la = int(255 * _ease_io(max(0.0, (t - 0.58) / 0.42)))
    _scene_text(draw, "Align",
                "Each copy is aligned to the template to identify amino-acid substitutions", la)
    return img


# ═══════════════════════════════════════════════════════════════════════════════
# SCENE 4 — REPORT
# ═══════════════════════════════════════════════════════════════════════════════

def scene_report(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)
    draw = ImageDraw.Draw(img)

    N      = len(MC_LABELS)
    CX     = 160
    CB     = 415
    CW     = 920
    CH     = 295
    CY     = CB - CH    # = 120
    MAX_V  = MC_COUNTS[0] * 1.12
    SLOT   = CW // N
    BAR_W  = SLOT - 18

    ax_a = int(255 * min(1.0, t * 5))

    # Axes
    draw.line([(CX, CY-4), (CX, CB+2)], fill=(*MID, ax_a), width=2)
    draw.line([(CX-2, CB), (CX+CW, CB)], fill=(*MID, ax_a), width=2)

    # Y ticks
    for count, lbl in [(50,"50"), (100,"100"), (150,"150"), (200,"200")]:
        ty = int(CB - count / MAX_V * CH)
        draw.line([(CX-5, ty), (CX+3, ty)], fill=(*MID, ax_a), width=1)
        tw = _textw(F_SM, lbl)
        draw.text((CX - tw - 7, ty - 8), lbl, font=F_SM, fill=(*MID, ax_a))

    draw.text((CX - 10, CY - 22), "read count", font=F_SM, fill=(*MID, ax_a))
    xlab = "amino-acid substitution"
    draw.text((CX + CW//2 - _textw(F_SM,xlab)//2, CB + 22), xlab,
              font=F_SM, fill=(*MID, ax_a))

    # Threshold line
    thresh_y = int(CB - MC_THRESH / MAX_V * CH)
    if ax_a > 60:
        draw.line([(CX, thresh_y), (CX+CW, thresh_y)],
                  fill=(*VIOLET, ax_a//2), width=1)
        lbl = f"threshold = {MC_THRESH}"
        draw.text((CX+CW - _textw(F_SM,lbl) - 4, thresh_y - 18),
                  lbl, font=F_SM, fill=(*VIOLET, ax_a))

    # Bars
    for k, (lbl, count) in enumerate(zip(MC_LABELS, MC_COUNTS)):
        delay = 0.04 * k
        bar_t = _ease_out(max(0.0, (t - delay - 0.04) / 0.96))
        bx    = CX + k * SLOT + 9
        bar_h = int(CH * count / MAX_V * bar_t)
        col   = CORAL if count >= MC_THRESH else SLATE

        if bar_h > 0:
            draw.rounded_rectangle((bx, CB-bar_h, bx+BAR_W, CB),
                                   radius=5, fill=(*col, 215))

        if bar_t > 0.6:
            ta  = int(215 * (bar_t - 0.6) / 0.4)
            num = str(count)
            tw  = _textw(F_SM, num)
            draw.text((bx + BAR_W//2 - tw//2, CB-bar_h-20),
                      num, font=F_SM, fill=(*INK, ta))

        # X-axis label (angled text not possible in PIL, draw straight)
        draw.text((bx + BAR_W//2 - _textw(F_XS,lbl)//2, CB+6),
                  lbl, font=F_XS, fill=(*MID, ax_a))

    # Legend
    if t > 0.70:
        la2 = int(220 * (t - 0.70) / 0.30)
        for col, txt, dy in [(CORAL,"≥ threshold (frequent)", 0),
                              (SLATE,"below threshold", 22)]:
            draw.rounded_rectangle((CX+8, CY+8+dy, CX+20, CY+20+dy),
                                   radius=3, fill=(*col, la2))
            draw.text((CX+26, CY+8+dy), txt, font=F_SM, fill=(*INK, la2))

    la = int(255 * _ease_io(max(0.0, (t - 0.58) / 0.42)))
    _scene_text(draw, "Report",
                "Substitution frequencies ranked — frequent variants highlighted above threshold", la)
    return img


# ═══════════════════════════════════════════════════════════════════════════════
# ASSEMBLE
# ═══════════════════════════════════════════════════════════════════════════════

def _crossfade(a: Image.Image, b: Image.Image, t: float) -> Image.Image:
    return Image.blend(a.convert("RGBA"), b.convert("RGBA"), _ease_io(t))


def build_frames() -> List[Image.Image]:
    fns      = [scene_locate, scene_extract, scene_align, scene_report]
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


_STATIC_COPY = Path("src/uht_tooling/web/static/animations/mutation_caller.gif")


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
