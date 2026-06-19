"""Profile Inserts — explainer GIF.

Scenes:
  1. Probe   — reads drop in; upstream/downstream probes glow amber; insert region labelled
  2. Extract — insert strips lift out of each read and stack below
  3. Profile — insert-length histogram grows; mean line and summary stats appear
  4. Report  — output files listed with key QC numbers
"""

from __future__ import annotations
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from PIL import Image, ImageDraw, ImageFilter, ImageFont

WIDTH, HEIGHT = 1200, 675
FPS           = 24
SCENE_FRAMES  = 48
FADE_FRAMES   = 12
HOLD_FRAMES   = 6

OUTPUT_PATH = Path("assets/animations/profile_inserts.gif")
POSTER_PATH = Path("assets/animations/profile_inserts_poster.png")

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
# [vl | upstream probe | INSERT | downstream probe | vr]
_VL, _UP, _INS, _DN, _VR = 0.12, 0.08, 0.56, 0.08, 0.16

UPS_START = _VL
UPS_END   = _VL + _UP
INS_START = UPS_END
INS_END   = UPS_END + _INS
DNS_START = INS_END
DNS_END   = INS_END + _DN

# ── Scene 3: length histogram data (simulated, bp) ────────────────────────────
# 20 bins, 50 bp wide, starting at 500 bp
HIST_BIN_START = 500
HIST_BIN_STEP  = 50
HIST_HEIGHTS   = [2, 4, 10, 22, 48, 88, 148, 210, 258, 280, 265, 222, 170, 118, 72, 40, 20, 8, 3, 1]
HIST_MEAN_BIN  = 9.0   # fractional bin index corresponding to mean (~950 bp)
HIST_MED_BIN   = 8.8

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
        "vl":  (0.0,       UPS_START),
        "up":  (UPS_START, UPS_END),
        "ins": (INS_START, INS_END),
        "dn":  (DNS_START, DNS_END),
        "vr":  (DNS_END,   1.0),
    }
    return {
        name: (x + 4 + int(inner * s), x + 4 + int(inner * e))
        for name, (s, e) in segs.items()
    }


def _draw_read_bar(
    img: Image.Image, x: int, y: int, w: int, h: int,
    alpha: int = 255,
    probe_glow: float = 0.0,
    ins_highlight: bool = False,
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

    seg("vl",  SLATE)
    seg("vr",  SLATE)
    seg("ins", GREEN)
    seg("up",  AMBER)
    seg("dn",  AMBER)

    # Insert label
    ix0, ix1 = sx["ins"]
    lbl = "insert"
    tw  = _textw(F_TINY, lbl)
    if ix1 - ix0 > tw + 4:
        draw.text((ix0 + (ix1-ix0)//2 - tw//2, y + h//2 - 6),
                  lbl, font=F_TINY, fill=(255, 255, 255, alpha))

    # Probe glow (both ends)
    if probe_glow > 0.01:
        ga = int(80 * probe_glow)
        for sn in ("up", "dn"):
            fx0, fx1 = sx[sn]
            ov = Image.new("RGBA", img.size, (0, 0, 0, 0))
            ImageDraw.Draw(ov).rounded_rectangle(
                (fx0-5, y-5, fx1+5, y+h+5), radius=r+5, fill=(*AMBER, ga))
            img.alpha_composite(ov.filter(ImageFilter.GaussianBlur(7)))

    # Insert highlight border
    if ins_highlight:
        ix0, ix1 = sx["ins"]
        ov = Image.new("RGBA", img.size, (0, 0, 0, 0))
        ImageDraw.Draw(ov).rounded_rectangle(
            (ix0-3, y+1, ix1+3, y+h-1), radius=r+3, fill=(*BLUE, 18))
        img.alpha_composite(ov.filter(ImageFilter.GaussianBlur(5)))
        ImageDraw.Draw(img).rounded_rectangle(
            (ix0-3, y+1, ix1+3, y+h-1), radius=r+3,
            outline=(*BLUE, min(255, alpha)), width=2)


def _draw_insert_strip(
    img: Image.Image, x: int, y: int, w: int, h: int,
    alpha: int = 255,
) -> None:
    draw = ImageDraw.Draw(img)
    draw.rounded_rectangle((x, y, x+w, y+h), radius=h//2,
                            fill=(*GREEN, alpha), outline=(*BORD, alpha//2), width=1)


# ══════════════════════════════════════════════════════════════════════════════
# SCENE 1 — PROBE
# ══════════════════════════════════════════════════════════════════════════════

def scene_probe(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)

    NUM  = 6
    RX, RW, RH, RSP = 90, 1010, 30, 56

    probe_glow = _ease_out(max(0.0, (t - 0.42) / 0.40))
    ins_hl     = t > 0.62

    for idx in range(NUM):
        delay   = idx * 0.07
        local_t = max(0.0, min(1.0, (t - delay) / max(1e-4, 1.0 - delay)))
        anim_t  = _elastic_out(local_t)
        y_final = 55 + idx * RSP
        y_now   = int(y_final - 80 * (1.0 - anim_t))
        a       = int(255 * min(1.0, local_t * 4))
        _draw_read_bar(img, RX, y_now, RW, RH, alpha=a,
                       probe_glow=probe_glow * min(1.0, local_t * 2),
                       ins_highlight=ins_hl)

    draw = ImageDraw.Draw(img)
    sx = _seg_xs(RX, RW)
    y0 = 55

    # Probe labels above first read
    if t > 0.52:
        la = int(220 * min(1.0, (t - 0.52) / 0.32))

        # Upstream probe label
        ux0, ux1 = sx["up"]
        ucx = (ux0 + ux1) // 2
        draw.line([(ucx, y0 - 4), (ucx, y0 - 20)], fill=(*AMBER, la), width=2)
        lbl = "upstream probe"
        draw.text((ucx - _textw(F_SM, lbl)//2, y0 - 38),
                  lbl, font=F_SM, fill=(*AMBER, la))

        # Downstream probe label
        dx0, dx1 = sx["dn"]
        dcx = (dx0 + dx1) // 2
        draw.line([(dcx, y0 - 4), (dcx, y0 - 20)], fill=(*AMBER, la), width=2)
        lbl2 = "downstream probe"
        draw.text((dcx - _textw(F_SM, lbl2)//2, y0 - 38),
                  lbl2, font=F_SM, fill=(*AMBER, la))

    # Insert bracket
    if t > 0.65:
        la = int(220 * min(1.0, (t - 0.65) / 0.28))
        ix0, ix1 = sx["ins"]
        gcx = (ix0 + ix1) // 2
        draw.line([(gcx, y0 - 4), (gcx, y0 - 20)], fill=(*GREEN, la), width=2)
        lbl = "insert sequence"
        draw.text((gcx - _textw(F_SM, lbl)//2, y0 - 38),
                  lbl, font=F_SM, fill=(*GREEN, la))
        if ins_hl:
            draw.line([(ix0, y0-4), (ix0, y0-16)], fill=(*BLUE, la), width=2)
            draw.line([(ix1, y0-4), (ix1, y0-16)], fill=(*BLUE, la), width=2)
            draw.line([(ix0, y0-16), (ix1, y0-16)], fill=(*BLUE, la), width=2)

    la = int(255 * _ease_io(max(0.0, (t - 0.58) / 0.42)))
    _scene_text(draw, "Probe",
                "Upstream and downstream probe pairs define the insert boundary in every read", la)
    return img


# ══════════════════════════════════════════════════════════════════════════════
# SCENE 2 — EXTRACT
# ══════════════════════════════════════════════════════════════════════════════

def scene_extract(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)

    NUM    = 6
    RX, RW, RH, RSP = 90, 1010, 22, 44

    # Static dimmed reads at top
    for idx in range(NUM):
        y = 42 + idx * RSP
        a = int(130 * (1.0 - min(1.0, t * 1.5)))
        if a > 5:
            _draw_read_bar(img, RX, y, RW, RH, alpha=a)

    # Insert strips fly out and land in stack
    sx_map  = _seg_xs(RX, RW)
    ix0_src, ix1_src = sx_map["ins"]

    STACK_X  = 220
    STACK_W  = RW - 340
    STACK_H  = 26
    STACK_SP = 32

    draw = ImageDraw.Draw(img)

    for idx in range(NUM):
        delay   = idx * 0.07
        local_t = _ease_out(max(0.0, (t - delay - 0.05) / max(1e-4, 0.95 - delay)))
        a       = int(255 * min(1.0, local_t * 3))
        if a < 5:
            continue

        src_y  = 42 + idx * RSP + (RH - STACK_H) // 2
        dst_y  = 240 + idx * STACK_SP
        y_cur  = int(src_y + (dst_y - src_y) * local_t)

        src_x  = ix0_src
        dst_x  = STACK_X
        x_cur  = int(src_x + (dst_x - src_x) * local_t)
        w_cur  = STACK_W if local_t > 0.8 else int(ix1_src - ix0_src)

        _draw_insert_strip(img, x_cur, y_cur, w_cur, STACK_H, alpha=a)

    if t > 0.55:
        la  = int(220 * (t - 0.55) / 0.45)
        draw = ImageDraw.Draw(img)
        lbl = "extracted inserts"
        draw.text((STACK_X + STACK_W + 16, 240 + NUM * STACK_SP // 2 - 10),
                  lbl, font=F_MED, fill=(*GREEN, la))

    draw = ImageDraw.Draw(img)
    la = int(255 * _ease_io(max(0.0, (t - 0.58) / 0.42)))
    _scene_text(draw, "Extract",
                "Insert sequences are pulled from every read spanning both probes", la)
    return img


# ══════════════════════════════════════════════════════════════════════════════
# SCENE 3 — PROFILE
# ══════════════════════════════════════════════════════════════════════════════

def scene_profile(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)
    draw = ImageDraw.Draw(img)

    N      = len(HIST_HEIGHTS)
    CX     = 140
    CB     = 420
    CW     = 740
    CH     = 290
    CY     = CB - CH
    MAX_H  = max(HIST_HEIGHTS) * 1.08
    SLOT   = CW // N
    BAR_W  = SLOT - 6

    ax_a = int(255 * min(1.0, t * 5))

    # Axes
    draw.line([(CX, CY-4),  (CX, CB+2)],   fill=(*MID, ax_a), width=2)
    draw.line([(CX-2, CB),  (CX+CW, CB)],  fill=(*MID, ax_a), width=2)

    # Y-axis ticks
    for count, lbl in [(50,"50"), (100,"100"), (150,"150"), (200,"200"), (250,"250")]:
        ty = int(CB - count / MAX_H * CH)
        draw.line([(CX-5, ty), (CX+3, ty)], fill=(*MID, ax_a), width=1)
        tw = _textw(F_SM, lbl)
        draw.text((CX - tw - 7, ty - 8), lbl, font=F_SM, fill=(*MID, ax_a))

    draw.text((CX - 10, CY - 22), "frequency", font=F_SM, fill=(*MID, ax_a))

    # X-axis label
    xlab = "insert length (bp)"
    draw.text((CX + CW//2 - _textw(F_SM, xlab)//2, CB + 22),
              xlab, font=F_SM, fill=(*MID, ax_a))

    # X-axis tick labels (every 5 bins = 250 bp)
    for i in range(0, N+1, 4):
        bx   = CX + i * SLOT
        lbl  = str(HIST_BIN_START + i * HIST_BIN_STEP)
        draw.text((bx - _textw(F_XS, lbl)//2, CB + 6), lbl, font=F_XS, fill=(*MID, ax_a))

    # Bars
    for k, h in enumerate(HIST_HEIGHTS):
        delay = 0.03 * k
        bar_t = _ease_out(max(0.0, (t - delay - 0.03) / 0.97))
        bx    = CX + k * SLOT + 3
        bar_h = int(CH * h / MAX_H * bar_t)
        if bar_h > 0:
            draw.rounded_rectangle((bx, CB-bar_h, bx+BAR_W, CB),
                                   radius=4, fill=(*BLUE, 210))

    # Mean line
    if t > 0.50:
        ma  = int(220 * min(1.0, (t - 0.50) / 0.30))
        mean_x = int(CX + HIST_MEAN_BIN * SLOT)
        draw.line([(mean_x, CY), (mean_x, CB)],
                  fill=(*CORAL, ma), width=2)
        # dashes drawn as line with short segments
        for dy in range(CY, CB, 10):
            draw.line([(mean_x, dy), (mean_x, min(dy+6, CB))],
                      fill=(*CORAL, ma), width=2)
        lbl = "mean ~950 bp"
        draw.text((mean_x + 6, CY + 4), lbl, font=F_XS, fill=(*CORAL, ma))

    # Stats panel (right side)
    if t > 0.62:
        sa  = int(210 * min(1.0, (t - 0.62) / 0.28))
        px  = CX + CW + 36
        py  = CY + 10
        stats = [
            ("Total inserts",    "4,712"),
            ("Mean length",      "950 bp"),
            ("Median length",    "940 bp"),
            ("GC content",       "52.3 %"),
            ("Unique sequences", "4,489"),
            ("Duplicate rate",   "4.7 %"),
        ]
        draw.text((px, py - 20), "QC Summary", font=F_SM, fill=(*MID, sa))
        for i, (key, val) in enumerate(stats):
            ky = py + i * 38
            draw.text((px, ky),      key, font=F_XS, fill=(*MID, sa))
            draw.text((px, ky + 14), val, font=F_SM, fill=(*INK, sa))
            draw.line([(px, ky + 34), (px + 220, ky + 34)],
                      fill=(*BORD, sa//2), width=1)

    la = int(255 * _ease_io(max(0.0, (t - 0.58) / 0.42)))
    _scene_text(draw, "Profile",
                "Insert-length distribution and QC metrics are computed across all extracted reads", la)
    return img


# ══════════════════════════════════════════════════════════════════════════════
# SCENE 4 — REPORT
# ══════════════════════════════════════════════════════════════════════════════

def scene_report(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)
    draw = ImageDraw.Draw(img)

    # Three output file cards arranged horizontally
    outputs = [
        {
            "file":  "extracted_inserts.fasta",
            "desc":  "All recovered insert sequences",
            "col":   GREEN,
            "stats": [("Inserts", "4,712"), ("Unique", "4,489")],
        },
        {
            "file":  "qc_report.txt",
            "desc":  "Summary statistics & probe QC",
            "col":   VIOLET,
            "stats": [("Mean len", "950 bp"), ("GC", "52.3 %")],
        },
        {
            "file":  "qc_plots.png",
            "desc":  "Multi-panel QC visualization",
            "col":   AMBER,
            "stats": [("12 panels", ""), ("probe perf", "")],
        },
    ]

    CARD_W = 300
    CARD_H = 240
    GAP    = 36
    TOTAL  = len(outputs) * CARD_W + (len(outputs) - 1) * GAP
    CX0    = (WIDTH - TOTAL) // 2
    CY0    = 90

    for i, out in enumerate(outputs):
        delay = 0.12 * i
        ct    = _ease_out(max(0.0, (t - delay) / max(1e-4, 1.0 - delay)))
        ca    = int(255 * min(1.0, ct * 4))
        if ca < 4:
            continue

        cx = CX0 + i * (CARD_W + GAP)
        cy = int(CY0 + 30 * (1.0 - ct))

        # Card shadow + border
        sh = Image.new("RGBA", img.size, (0, 0, 0, 0))
        sd = ImageDraw.Draw(sh)
        sd.rounded_rectangle((cx+3, cy+6, cx+CARD_W+3, cy+CARD_H+6),
                              radius=14, fill=(0, 0, 0, 24))
        img.alpha_composite(sh.filter(ImageFilter.GaussianBlur(8)))
        ImageDraw.Draw(img).rounded_rectangle(
            (cx, cy, cx+CARD_W, cy+CARD_H), radius=14,
            fill=(*CARD, ca), outline=(*BORD, ca), width=1)

        draw = ImageDraw.Draw(img)

        # Colour strip at top
        draw.rounded_rectangle((cx, cy, cx+CARD_W, cy+8),
                                radius=14, fill=(*out["col"], ca))
        draw.rectangle((cx, cy+4, cx+CARD_W, cy+8), fill=(*out["col"], ca))

        # Filename
        draw.text((cx + 16, cy + 24), out["file"],
                  font=F_SM, fill=(*INK, ca))

        # Description
        desc_y = cy + 50
        draw.text((cx + 16, desc_y), out["desc"],
                  font=F_XS, fill=(*MID, ca))

        # Divider
        draw.line([(cx+16, desc_y+22), (cx+CARD_W-16, desc_y+22)],
                  fill=(*BORD, ca//2), width=1)

        # Stats rows
        for j, (key, val) in enumerate(out["stats"]):
            sy = desc_y + 34 + j * 44
            draw.text((cx + 16, sy),      key, font=F_XS, fill=(*MID, ca))
            if val:
                draw.text((cx + 16, sy + 14), val, font=F_SM, fill=(*out["col"], ca))

        # Large file-type icon (letter in circle)
        icon_ext = out["file"].rsplit(".", 1)[-1].upper()
        ico_cx   = cx + CARD_W - 44
        ico_cy   = cy + CARD_H - 44
        draw.ellipse((ico_cx-18, ico_cy-18, ico_cx+18, ico_cy+18),
                     fill=(*out["col"], ca//3))
        iw = _textw(F_XS, icon_ext)
        draw.text((ico_cx - iw//2, ico_cy - 6), icon_ext,
                  font=F_XS, fill=(*out["col"], ca))

    la = int(255 * _ease_io(max(0.0, (t - 0.58) / 0.42)))
    _scene_text(draw, "Report",
                "Insert FASTA, QC statistics, and multi-panel plots written per sample", la)
    return img


# ══════════════════════════════════════════════════════════════════════════════
# ASSEMBLE
# ══════════════════════════════════════════════════════════════════════════════

def _crossfade(a: Image.Image, b: Image.Image, t: float) -> Image.Image:
    return Image.blend(a.convert("RGBA"), b.convert("RGBA"), _ease_io(t))


def build_frames() -> List[Image.Image]:
    fns      = [scene_probe, scene_extract, scene_profile, scene_report]
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


_STATIC_COPY = Path("src/uht_tooling/web/static/animations/profile_inserts.gif")


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
