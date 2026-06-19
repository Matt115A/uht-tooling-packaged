"""UMI Hunter — explainer GIF.

Scenes:
  1. Locate  — reads drop in; UMI flanks glow amber, gene flanks glow teal; both regions labelled
  2. Cluster — reads with similar UMI barcodes animate into two cluster groups
  3. Align   — cluster's gene copies align to the template; consensus sequence emerges
  4. Report  — cluster-size bar chart; output file listing fades in
"""

from __future__ import annotations
import math
from pathlib import Path
from typing import Dict, List, Tuple

from PIL import Image, ImageDraw, ImageFilter, ImageFont

WIDTH, HEIGHT = 1200, 675
FPS           = 24
SCENE_FRAMES  = 48
FADE_FRAMES   = 12
HOLD_FRAMES   = 6

OUTPUT_PATH = Path("assets/animations/umi_hunter.gif")
POSTER_PATH = Path("assets/animations/umi_hunter_poster.png")

# ── Colours ────────────────────────────────────────────────────────────────────
BG     = (245, 245, 247)
CARD   = (255, 255, 255)
INK    = (29,  29,  31)
MID    = (110, 110, 115)
BLUE   = (0,   113, 227)
BORD   = (210, 210, 215)
TEAL   = (0,   196, 175)
GREEN  = (40,  185, 115)
AMBER  = (255, 180,  50)
CORAL  = (255,  80,  65)
VIOLET = (155, 100, 240)
SLATE  = (75,  115, 160)

# ── Read anatomy fractions ─────────────────────────────────────────────────────
# [vl | ufs | umi | ufe | mid | gfs | gene | gfe | vr]
_VL, _UFS, _UMI, _UFE = 0.07, 0.04, 0.08, 0.04
_MID, _GFS, _GENE, _GFE, _VR = 0.10, 0.04, 0.47, 0.04, 0.12

UFS_START  = _VL
UMI_START  = _VL + _UFS
UMI_END    = UMI_START + _UMI
UFE_END    = UMI_END + _UFE
GFS_START  = UFE_END + _MID
GFS_END    = GFS_START + _GFS
GENE_START = GFS_END
GENE_END   = GENE_START + _GENE
GFE_END    = GENE_END + _GFE

# ── Cluster membership (5 reads: 0→A, 1→A, 2→B, 3→A, 4→B) ────────────────────
READ_CLUSTERS = [0, 0, 1, 0, 1]
CLUSTER_COLS  = [BLUE, CORAL]
CLUSTER_NAMES = ["Cluster A  (n = 3)", "Cluster B  (n = 2)"]

# ── Report scene ───────────────────────────────────────────────────────────────
CLUST_LABELS = ["Clust_1", "Clust_2", "Clust_3", "Clust_4", "Clust_5", "Clust_6"]
CLUST_COUNTS = [312, 187, 143, 95, 61, 28]
MIN_CLUSTER  = 50

_FONT_PATH = "/System/Library/Fonts/Supplemental/Arial.ttf"
_BOLD_PATH = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"


def _font(size: int, bold: bool = False):
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
F_XS   = _font(10)


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


def _base_img() -> Image.Image:
    return Image.new("RGBA", (WIDTH, HEIGHT), BG)

def _draw_card(img: Image.Image, x0: int, y0: int, x1: int, y1: int) -> None:
    sh = Image.new("RGBA", img.size, (0, 0, 0, 0))
    sd = ImageDraw.Draw(sh)
    sd.rounded_rectangle((x0+4, y0+8, x1+4, y1+8), radius=18, fill=(0, 0, 0, 32))
    img.alpha_composite(sh.filter(ImageFilter.GaussianBlur(12)))
    ImageDraw.Draw(img).rounded_rectangle(
        (x0, y0, x1, y1), radius=18,
        fill=(*CARD, 255), outline=(*BORD, 255), width=1)

def _scene_text(draw: ImageDraw.ImageDraw, title: str, subtitle: str, alpha: int) -> None:
    draw.rectangle((70, 488, WIDTH-70, 490), fill=(*BLUE, alpha//3))
    draw.text((80, 502),  title,    font=F_HUGE, fill=(*INK, alpha))
    draw.text((82, 588),  subtitle, font=F_MED,  fill=(*MID, alpha))

def _textw(font, s: str) -> int:
    try:
        return int(font.getlength(s))
    except AttributeError:
        return len(s) * 8


def _seg_xs(x: int, w: int) -> Dict[str, Tuple[int, int]]:
    inner = w - 8
    segs = {
        "vl":   (0.0,        UFS_START),
        "ufs":  (UFS_START,  UMI_START),
        "umi":  (UMI_START,  UMI_END),
        "ufe":  (UMI_END,    UFE_END),
        "mid":  (UFE_END,    GFS_START),
        "gfs":  (GFS_START,  GFS_END),
        "gene": (GENE_START, GENE_END),
        "gfe":  (GENE_END,   GFE_END),
        "vr":   (GFE_END,    1.0),
    }
    return {
        name: (x + 4 + int(inner * s), x + 4 + int(inner * e))
        for name, (s, e) in segs.items()
    }


def _draw_read_bar(
    img: Image.Image, x: int, y: int, w: int, h: int,
    alpha: int = 255,
    umi_glow: float = 0.0,
    gene_glow: float = 0.0,
    umi_highlight: bool = False,
    gene_highlight: bool = False,
    umi_col: Tuple[int, int, int] = VIOLET,
) -> None:
    draw = ImageDraw.Draw(img)
    draw.rounded_rectangle((x, y, x+w, y+h), radius=h//2,
                            fill=(215, 222, 232, alpha), outline=(*BORD, alpha), width=1)
    sx = _seg_xs(x, w)
    inner_h = h - 8
    r = max(2, inner_h // 2)

    def seg(name: str, col: Tuple[int, int, int]) -> None:
        x0, x1 = sx[name]
        if x1 > x0 + 1:
            draw.rounded_rectangle((x0, y+4, x1, y+h-4), radius=r, fill=(*col, alpha))

    seg("vl",   SLATE)
    seg("mid",  SLATE)
    seg("vr",   SLATE)
    seg("gene", GREEN)
    seg("gfs",  TEAL)
    seg("gfe",  TEAL)
    seg("umi",  umi_col)
    seg("ufs",  AMBER)
    seg("ufe",  AMBER)

    # Gene label
    gx0, gx1 = sx["gene"]
    lbl = "gene"
    tw  = _textw(F_TINY, lbl)
    if gx1 - gx0 > tw + 4:
        draw.text((gx0 + (gx1-gx0)//2 - tw//2, y + h//2 - 6),
                  lbl, font=F_TINY, fill=(255, 255, 255, alpha))

    # UMI label
    ux0, ux1 = sx["umi"]
    lbl2 = "UMI"
    tw2  = _textw(F_TINY, lbl2)
    if ux1 - ux0 > tw2 + 2:
        draw.text((ux0 + (ux1-ux0)//2 - tw2//2, y + h//2 - 6),
                  lbl2, font=F_TINY, fill=(255, 255, 255, alpha))

    # UMI glow
    if umi_glow > 0.01:
        ga = int(80 * umi_glow)
        for sn, col_g in [("ufs", AMBER), ("ufe", AMBER), ("umi", umi_col)]:
            fx0, fx1 = sx[sn]
            ov = Image.new("RGBA", img.size, (0, 0, 0, 0))
            ImageDraw.Draw(ov).rounded_rectangle(
                (fx0-5, y-5, fx1+5, y+h+5), radius=r+5, fill=(*col_g, ga))
            img.alpha_composite(ov.filter(ImageFilter.GaussianBlur(7)))

    # Gene glow
    if gene_glow > 0.01:
        ga = int(80 * gene_glow)
        for sn in ("gfs", "gfe"):
            fx0, fx1 = sx[sn]
            ov = Image.new("RGBA", img.size, (0, 0, 0, 0))
            ImageDraw.Draw(ov).rounded_rectangle(
                (fx0-5, y-5, fx1+5, y+h+5), radius=r+5, fill=(*TEAL, ga))
            img.alpha_composite(ov.filter(ImageFilter.GaussianBlur(7)))

    # UMI highlight border
    if umi_highlight:
        ov = Image.new("RGBA", img.size, (0, 0, 0, 0))
        ImageDraw.Draw(ov).rounded_rectangle(
            (ux0-3, y+1, ux1+3, y+h-1), radius=r+3, fill=(*umi_col, 20))
        img.alpha_composite(ov.filter(ImageFilter.GaussianBlur(5)))
        ImageDraw.Draw(img).rounded_rectangle(
            (ux0-3, y+1, ux1+3, y+h-1), radius=r+3,
            outline=(*umi_col, min(255, alpha)), width=2)

    # Gene highlight border
    if gene_highlight:
        ov = Image.new("RGBA", img.size, (0, 0, 0, 0))
        ImageDraw.Draw(ov).rounded_rectangle(
            (gx0-3, y+1, gx1+3, y+h-1), radius=r+3, fill=(*GREEN, 20))
        img.alpha_composite(ov.filter(ImageFilter.GaussianBlur(5)))
        ImageDraw.Draw(img).rounded_rectangle(
            (gx0-3, y+1, gx1+3, y+h-1), radius=r+3,
            outline=(*GREEN, min(255, alpha)), width=2)


def _draw_gene_strip(
    img: Image.Image, x: int, y: int, w: int, h: int,
    alpha: int = 255,
    col: Tuple[int, int, int] = GREEN,
) -> None:
    draw = ImageDraw.Draw(img)
    draw.rounded_rectangle((x, y, x+w, y+h), radius=h//2,
                            fill=(*col, alpha), outline=(*BORD, alpha//2), width=1)


# ══════════════════════════════════════════════════════════════════════════════
# SCENE 1 — LOCATE
# ══════════════════════════════════════════════════════════════════════════════

def scene_locate(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)

    NUM  = 5
    RX, RW, RH, RSP = 90, 1010, 30, 56

    umi_glow  = _ease_out(max(0.0, (t - 0.42) / 0.40))
    gene_glow = _ease_out(max(0.0, (t - 0.54) / 0.35))
    umi_hl    = t > 0.60
    gene_hl   = t > 0.70

    for idx in range(NUM):
        delay   = idx * 0.08
        local_t = max(0.0, min(1.0, (t - delay) / max(1e-4, 1.0 - delay)))
        anim_t  = _elastic_out(local_t)
        y_final = 62 + idx * RSP
        y_now   = int(y_final - 80 * (1.0 - anim_t))
        a       = int(255 * min(1.0, local_t * 4))
        _draw_read_bar(
            img, RX, y_now, RW, RH, alpha=a,
            umi_glow=umi_glow * min(1.0, local_t * 2),
            gene_glow=gene_glow * min(1.0, local_t * 2),
            umi_highlight=umi_hl,
            gene_highlight=gene_hl,
        )

    draw = ImageDraw.Draw(img)
    sx = _seg_xs(RX, RW)
    y0 = 62  # first read top y

    # UMI region annotation
    if t > 0.52:
        la  = int(220 * min(1.0, (t - 0.52) / 0.32))
        ux0, ux1 = sx["umi"]
        ucx = (ux0 + ux1) // 2
        draw.line([(ucx, y0 - 4), (ucx, y0 - 20)], fill=(*VIOLET, la), width=2)
        lbl = "UMI barcode"
        draw.text((ucx - _textw(F_SM, lbl)//2, y0 - 38),
                  lbl, font=F_SM, fill=(*VIOLET, la))
        if umi_hl:
            draw.line([(ux0, y0-4), (ux0, y0-16)], fill=(*VIOLET, la), width=2)
            draw.line([(ux1, y0-4), (ux1, y0-16)], fill=(*VIOLET, la), width=2)
            draw.line([(ux0, y0-16), (ux1, y0-16)], fill=(*VIOLET, la), width=2)

    # Gene region annotation
    if t > 0.63:
        la  = int(220 * min(1.0, (t - 0.63) / 0.28))
        gx0, gx1 = sx["gene"]
        gcx = (gx0 + gx1) // 2
        draw.line([(gcx, y0 - 4), (gcx, y0 - 20)], fill=(*GREEN, la), width=2)
        lbl = "gene insert"
        draw.text((gcx - _textw(F_SM, lbl)//2, y0 - 38),
                  lbl, font=F_SM, fill=(*GREEN, la))
        if gene_hl:
            draw.line([(gx0, y0-4), (gx0, y0-16)], fill=(*BLUE, la), width=2)
            draw.line([(gx1, y0-4), (gx1, y0-16)], fill=(*BLUE, la), width=2)
            draw.line([(gx0, y0-16), (gx1, y0-16)], fill=(*BLUE, la), width=2)

    la = int(255 * _ease_io(max(0.0, (t - 0.58) / 0.42)))
    _scene_text(draw, "Locate",
                "Flanking sequences pinpoint the UMI barcode and gene insert in every read", la)
    return img


# ══════════════════════════════════════════════════════════════════════════════
# SCENE 2 — CLUSTER
# ══════════════════════════════════════════════════════════════════════════════

def scene_cluster(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)

    NUM  = 5
    RH, RSP = 26, 50
    RW_FULL = 980
    RX_FULL = 110  # centered x when stacked

    # Left group (cluster A) and right group (cluster B) final geometry
    LX, LW = 100, 460
    RX_GRP, RW_GRP = 620, 460

    # Within each group, how many reads
    a_indices = [i for i, c in enumerate(READ_CLUSTERS) if c == 0]  # [0,1,3]
    b_indices = [i for i, c in enumerate(READ_CLUSTERS) if c == 1]  # [2,4]

    t_split = 0.30  # when reads start separating

    draw = ImageDraw.Draw(img)

    for idx in range(NUM):
        cluster   = READ_CLUSTERS[idx]
        umi_col   = CLUSTER_COLS[cluster]
        a         = int(255 * min(1.0, t * 10))

        # Initial position (single centered column)
        y_init = 68 + idx * RSP

        # Final position (in their respective group column)
        if cluster == 0:
            gi      = a_indices.index(idx)
            x_fin   = LX
            w_fin   = LW
            y_fin   = 68 + gi * (RSP + 10)
        else:
            gi      = b_indices.index(idx)
            x_fin   = RX_GRP
            w_fin   = RW_GRP
            y_fin   = 68 + gi * (RSP + 10)

        sep_t = _ease_out(max(0.0, (t - t_split) / (1.0 - t_split)))

        x_now = int(RX_FULL + (x_fin - RX_FULL) * sep_t)
        w_now = int(RW_FULL + (w_fin - RW_FULL) * sep_t)
        y_now = int(y_init  + (y_fin  - y_init)  * sep_t)

        _draw_read_bar(
            img, x_now, y_now, w_now, RH, alpha=a,
            umi_highlight=True,
            umi_col=umi_col,
        )

    # Cluster labels (appear late)
    if t > 0.65:
        la = int(220 * min(1.0, (t - 0.65) / 0.30))

        # Cluster A label (left)
        a_last_y = 68 + (len(a_indices) - 1) * (RSP + 10)
        mid_y_a  = (68 + a_last_y + RH) // 2
        lbl_a    = CLUSTER_NAMES[0]
        tw_a     = _textw(F_SM, lbl_a)
        draw.text((LX + LW//2 - tw_a//2, a_last_y + RH + 12),
                  lbl_a, font=F_SM, fill=(*BLUE, la))
        # Bracket around cluster A
        y_top_a, y_bot_a = 60, a_last_y + RH + 8
        draw.rounded_rectangle(
            (LX - 8, y_top_a, LX + LW + 8, y_bot_a),
            radius=10, outline=(*BLUE, la//2), width=2)

        # Cluster B label (right)
        b_last_y = 68 + (len(b_indices) - 1) * (RSP + 10)
        lbl_b    = CLUSTER_NAMES[1]
        tw_b     = _textw(F_SM, lbl_b)
        draw.text((RX_GRP + RW_GRP//2 - tw_b//2, b_last_y + RH + 12),
                  lbl_b, font=F_SM, fill=(*CORAL, la))
        # Bracket around cluster B
        y_bot_b = b_last_y + RH + 8
        draw.rounded_rectangle(
            (RX_GRP - 8, y_top_a, RX_GRP + RW_GRP + 8, y_bot_b),
            radius=10, outline=(*CORAL, la//2), width=2)

    la = int(255 * _ease_io(max(0.0, (t - 0.58) / 0.42)))
    _scene_text(draw, "Cluster",
                "Reads sharing near-identical UMI barcodes are grouped into clusters", la)
    return img


# ══════════════════════════════════════════════════════════════════════════════
# SCENE 3 — ALIGN
# ══════════════════════════════════════════════════════════════════════════════

def scene_align(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)

    AX, AW, AH = 200, 800, 28
    NUM_READS   = 3
    RSP         = 50

    draw = ImageDraw.Draw(img)

    # Template bar at top
    tmpl_a = int(255 * min(1.0, t * 8))
    draw.rounded_rectangle((AX, 52, AX+AW, 52+AH), radius=AH//2,
                            fill=(*BLUE, tmpl_a), outline=(*BORD, tmpl_a//2), width=1)
    if tmpl_a > 60:
        lbl = "template"
        tw  = _textw(F_SM, lbl)
        draw.text((AX - tw - 10, 52 + AH//2 - 8), lbl,
                  font=F_SM, fill=(*BLUE, tmpl_a))

    # Gene strips from cluster A slide in and align
    for ri in range(NUM_READS):
        delay   = ri * 0.07 + 0.10
        local_t = _ease_out(max(0.0, (t - delay) / max(1e-4, 1.0 - delay)))
        a       = int(255 * min(1.0, local_t * 4))
        if a < 4:
            continue

        y_final = 52 + AH + 22 + ri * RSP
        x_start = AX + 300
        x_cur   = int(x_start + (AX - x_start) * local_t)

        _draw_gene_strip(img, x_cur, y_final, AW, AH, alpha=a, col=GREEN)

        if local_t > 0.5:
            la2 = int(a * (local_t - 0.5) / 0.5)
            draw = ImageDraw.Draw(img)
            lbl = f"read {ri + 1}"
            draw.text((x_cur - _textw(F_SM, lbl) - 10, y_final + AH//2 - 8),
                      lbl, font=F_SM, fill=(*MID, la2))

    # Consensus bar grows from left
    cons_y = 52 + AH + 22 + NUM_READS * RSP + 18
    cons_t = _ease_out(max(0.0, (t - 0.60) / 0.40))
    if cons_t > 0.01:
        ca = int(255 * min(1.0, cons_t * 3))
        draw = ImageDraw.Draw(img)
        cons_w = int(AW * cons_t)
        draw.rounded_rectangle((AX, cons_y, AX + cons_w, cons_y + AH),
                                radius=AH//2, fill=(*TEAL, ca),
                                outline=(*BORD, ca//2), width=1)
        if cons_t > 0.4:
            la3 = int(ca * (cons_t - 0.4) / 0.6)
            lbl = "consensus"
            tw  = _textw(F_SM, lbl)
            draw.text((AX - tw - 10, cons_y + AH//2 - 8),
                      lbl, font=F_SM, fill=(*TEAL, la3))

    draw = ImageDraw.Draw(img)
    la = int(255 * _ease_io(max(0.0, (t - 0.58) / 0.42)))
    _scene_text(draw, "Align",
                "Gene copies within each cluster align to the template — consensus is called per position", la)
    return img


# ══════════════════════════════════════════════════════════════════════════════
# SCENE 4 — REPORT
# ══════════════════════════════════════════════════════════════════════════════

def scene_report(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)
    draw = ImageDraw.Draw(img)

    N      = len(CLUST_LABELS)
    CX     = 160
    CB     = 415
    CW     = 700
    CH     = 295
    CY     = CB - CH
    MAX_V  = CLUST_COUNTS[0] * 1.12
    SLOT   = CW // N
    BAR_W  = SLOT - 16

    ax_a = int(255 * min(1.0, t * 5))

    # Axes
    draw.line([(CX, CY-4), (CX, CB+2)],    fill=(*MID, ax_a), width=2)
    draw.line([(CX-2, CB), (CX+CW, CB)],   fill=(*MID, ax_a), width=2)

    # Y-axis ticks
    for count, lbl in [(50,"50"), (100,"100"), (150,"150"), (200,"200"), (250,"250"), (300,"300")]:
        ty = int(CB - count / MAX_V * CH)
        draw.line([(CX-5, ty), (CX+3, ty)], fill=(*MID, ax_a), width=1)
        tw = _textw(F_SM, lbl)
        draw.text((CX - tw - 7, ty - 8), lbl, font=F_SM, fill=(*MID, ax_a))

    draw.text((CX - 10, CY - 22), "read count", font=F_SM, fill=(*MID, ax_a))
    xlab = "UMI cluster"
    draw.text((CX + CW//2 - _textw(F_SM, xlab)//2, CB + 22),
              xlab, font=F_SM, fill=(*MID, ax_a))

    # Min-cluster threshold line
    thresh_y = int(CB - MIN_CLUSTER / MAX_V * CH)
    if ax_a > 60:
        draw.line([(CX, thresh_y), (CX+CW, thresh_y)],
                  fill=(*VIOLET, ax_a//2), width=1)
        lbl = f"min cluster = {MIN_CLUSTER}"
        draw.text((CX + CW - _textw(F_SM, lbl) - 4, thresh_y - 18),
                  lbl, font=F_SM, fill=(*VIOLET, ax_a))

    # Bars
    for k, (lbl, count) in enumerate(zip(CLUST_LABELS, CLUST_COUNTS)):
        delay = 0.05 * k
        bar_t = _ease_out(max(0.0, (t - delay - 0.04) / 0.96))
        bx    = CX + k * SLOT + 8
        bar_h = int(CH * count / MAX_V * bar_t)
        col   = BLUE if count >= MIN_CLUSTER else SLATE

        if bar_h > 0:
            draw.rounded_rectangle((bx, CB-bar_h, bx+BAR_W, CB),
                                   radius=5, fill=(*col, 215))

        if bar_t > 0.6:
            ta  = int(215 * (bar_t - 0.6) / 0.4)
            num = str(count)
            tw  = _textw(F_SM, num)
            draw.text((bx + BAR_W//2 - tw//2, CB - bar_h - 20),
                      num, font=F_SM, fill=(*INK, ta))

        draw.text((bx + BAR_W//2 - _textw(F_XS, lbl)//2, CB+6),
                  lbl, font=F_XS, fill=(*MID, ax_a))

    # Output file listing (right panel)
    if t > 0.55:
        fa  = int(210 * min(1.0, (t - 0.55) / 0.35))
        px  = CX + CW + 40
        py  = CY + 20
        out_files = [
            ("*_UMI_clusters.csv",    "raw UMI cluster table",      VIOLET),
            ("*_gene_consensus.csv",  "per-cluster consensus gene",  GREEN),
            ("*_consensuses.fasta",   "consensus FASTA sequences",   TEAL),
        ]
        draw.text((px, py - 22), "Outputs per sample:", font=F_SM, fill=(*MID, fa))
        for i, (fname, desc, col) in enumerate(out_files):
            fy = py + i * 58
            draw.rounded_rectangle((px, fy, px + 14, fy + 14),
                                   radius=3, fill=(*col, fa))
            draw.text((px + 20, fy - 1), fname, font=F_SM, fill=(*INK, fa))
            draw.text((px + 20, fy + 18), desc, font=F_XS, fill=(*MID, fa))

    la = int(255 * _ease_io(max(0.0, (t - 0.58) / 0.42)))
    _scene_text(draw, "Report",
                "Cluster abundances, consensus genes, and FASTA sequences written per sample", la)
    return img


# ══════════════════════════════════════════════════════════════════════════════
# ASSEMBLE
# ══════════════════════════════════════════════════════════════════════════════

def _crossfade(a: Image.Image, b: Image.Image, t: float) -> Image.Image:
    return Image.blend(a.convert("RGBA"), b.convert("RGBA"), _ease_io(t))


def build_frames() -> List[Image.Image]:
    fns      = [scene_locate, scene_cluster, scene_align, scene_report]
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


_STATIC_COPY = Path("src/uht_tooling/web/static/animations/umi_hunter.gif")


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
