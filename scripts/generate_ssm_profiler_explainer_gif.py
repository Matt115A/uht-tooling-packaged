"""SSM Profiler — explainer GIF.

Scenes:
  1. Target   — ROI bar with target codon positions lighting up in gold
  2. Decode   — top amino acids per site drop in as coloured chips
  3. Heatmap  — full 20-AA × 3-site frequency heatmap fills in
  4. Assess   — per-read mutation-load distribution histogram
"""

from __future__ import annotations
import math
import random
from pathlib import Path
from typing import Dict, List, Tuple

from PIL import Image, ImageDraw, ImageFilter, ImageFont

WIDTH, HEIGHT = 1200, 675
FPS          = 24
SCENE_FRAMES = 48
FADE_FRAMES  = 12
HOLD_FRAMES  = 6

OUTPUT_PATH = Path("assets/animations/ssm_profiler.gif")
POSTER_PATH = Path("assets/animations/ssm_profiler_poster.png")

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
BAR_BG = (248, 248, 250)

# ── Construct segments (reused from synthetic gene pool style) ────────────────
SEG_FRACS  = [0.18, 0.28, 0.10, 0.07, 0.15, 0.22]
SEG_COLORS = [TEAL, GREEN, AMBER, CORAL, VIOLET, SLATE]

# ── Target-site positions as fraction within bar inner width ──────────────────
# Three target sites at 28 %, 52 %, 76 % of the inner bar extent
TARGET_SITES      = [45, 92, 163]
TARGET_FRACS      = [0.28, 0.52, 0.76]   # fraction across full inner bar
TARGET_SITE_NAMES = ["AA 45", "AA 92", "AA 163"]

# ── Amino-acid ordering and properties ───────────────────────────────────────
# Sorted by NNK codon count (high → low), then alphabetically
AA_ORDER = ["R", "S", "L", "A", "G", "P", "T", "V",
            "C", "D", "E", "F", "H", "I", "K", "M", "N", "Q", "W", "Y"]

# NNK codon counts (total 31, stop excluded)
NNK_COUNTS: Dict[str, int] = {
    "R": 3, "S": 3, "L": 3,
    "A": 2, "G": 2, "P": 2, "T": 2, "V": 2,
    "C": 1, "D": 1, "E": 1, "F": 1, "H": 1,
    "I": 1, "K": 1, "M": 1, "N": 1, "Q": 1, "W": 1, "Y": 1,
}
_NNK_TOTAL   = sum(NNK_COUNTS.values())  # 31
NNK_EXPECTED = {aa: n / _NNK_TOTAL for aa, n in NNK_COUNTS.items()}

# Per-site observed AA distributions (with realistic noise on NNK)
_rnd = random.Random(17)

def _make_obs(missing: List[str] = ()) -> Dict[str, float]:
    obs: Dict[str, float] = {}
    for aa in AA_ORDER:
        noise   = _rnd.gauss(0, 0.008)
        obs[aa] = max(0.0, NNK_EXPECTED[aa] + noise)
    for aa in missing:
        obs[aa] = max(0.0, obs[aa] * 0.08)
    total = sum(obs.values())
    return {aa: v / total for aa, v in obs.items()}

SSM_OBSERVED = [
    _make_obs(),                  # site 45  — good coverage
    _make_obs(missing=["C","W"]), # site 92  — missing two rare codons
    _make_obs(),                  # site 163 — good coverage
]

# Top-6 AAs per site (for Scene 2 chips)
TOP6 = [
    sorted(obs.items(), key=lambda x: -x[1])[:6]
    for obs in SSM_OBSERVED
]

# ── AA chemistry colours ──────────────────────────────────────────────────────
_AA_COL: Dict[str, Tuple[int,int,int]] = {
    "R": BLUE, "K": BLUE, "H": (0, 150, 200),
    "D": CORAL, "E": CORAL,
    "S": TEAL, "T": TEAL, "C": TEAL, "Y": TEAL, "N": TEAL, "Q": TEAL,
    "G": AMBER,
}
_HYDROPHOBIC = {"A","V","I","L","M","F","W","P"}
for _aa in _HYDROPHOBIC:
    _AA_COL[_aa] = SLATE

# ── Mutation-load distribution (scene 4) ─────────────────────────────────────
LOAD_VALS   = [0.08, 0.28, 0.37, 0.27]  # P(0), P(1), P(2), P(3) sites mutated
LOAD_LABELS = ["0", "1", "2", "3"]
LOAD_COLS   = [SLATE, AMBER, GREEN, TEAL]

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

def _draw_roi_bar(img: Image.Image, x: int, y: int, w: int, h: int,
                  alpha: int = 255, highlight_targets: bool = False,
                  target_t: float = 0.0) -> None:
    """Draw the gene construct bar (reusing synthetic-gene-pool visual language)."""
    draw = ImageDraw.Draw(img)
    draw.rounded_rectangle((x, y, x+w, y+h), radius=h//2,
                            fill=(*BAR_BG, alpha), outline=(*BORD, alpha), width=1)
    cursor, inner = x + 4, w - 8
    for i, (frac, col) in enumerate(zip(SEG_FRACS, SEG_COLORS)):
        sw = int(inner * frac)
        if i == len(SEG_FRACS) - 1:
            sw = (x + w - 4) - cursor
        if sw < 2:
            cursor += sw + 2
            continue
        draw.rounded_rectangle((cursor, y+4, cursor+sw, y+h-4),
                                radius=(h-8)//2, fill=(*col, alpha))
        if i == 1:  # gene segment label
            lbl = "gene"
            tw  = _textw(F_TINY, lbl)
            draw.text((cursor + sw//2 - tw//2, y+h//2-6), lbl,
                      font=F_TINY, fill=(255,255,255,min(255,alpha)))
        cursor += sw + 2

    if not highlight_targets or target_t < 0.01:
        return

    # Golden vertical bands at each target site
    ta = int(220 * target_t)
    for frac in TARGET_FRACS:
        tx   = x + 4 + int(inner * frac)
        band = int(6 + 8 * target_t)
        ov   = Image.new("RGBA", img.size, (0,0,0,0))
        ImageDraw.Draw(ov).rectangle(
            (tx - band, y - int(28*target_t),
             tx + band, y + h + int(28*target_t)),
            fill=(*AMBER, ta))
        img.alpha_composite(ov)


def _heatmap_color(freq: float, max_freq: float = 0.12) -> Tuple[int,int,int]:
    """White → amber → coral gradient for heatmap cells."""
    t = min(1.0, freq / max_freq)
    if t < 0.5:
        tt = t * 2
        return (
            int(255 + tt * (AMBER[0] - 255)),
            int(255 + tt * (AMBER[1] - 255)),
            int(255 + tt * (AMBER[2] - 255)),
        )
    tt = (t - 0.5) * 2
    return (
        int(AMBER[0] + tt * (CORAL[0] - AMBER[0])),
        int(AMBER[1] + tt * (CORAL[1] - AMBER[1])),
        int(AMBER[2] + tt * (CORAL[2] - AMBER[2])),
    )


# ═══════════════════════════════════════════════════════════════════════════════
# SCENE 1 — TARGET
# ═══════════════════════════════════════════════════════════════════════════════

def scene_target(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)

    BAR_X, BAR_Y, BAR_W, BAR_H = 110, 180, 980, 44

    tgt_t  = _ease_out(max(0.0, (t - 0.22) / 0.78))
    bar_a  = int(255 * min(1.0, t * 6))
    _draw_roi_bar(img, BAR_X, BAR_Y, BAR_W, BAR_H,
                  alpha=bar_a, highlight_targets=True, target_t=tgt_t)

    draw = ImageDraw.Draw(img)

    # Site labels below the bands
    if tgt_t > 0.3:
        la = int(220 * (tgt_t - 0.3) / 0.7)
        inner = BAR_W - 8
        for name, frac in zip(TARGET_SITE_NAMES, TARGET_FRACS):
            tx = BAR_X + 4 + int(inner * frac)
            tw = _textw(F_SM, name)
            draw.text((tx - tw//2, BAR_Y + BAR_H + 36), name,
                      font=F_SM, fill=(*AMBER, la))
            draw.line([(tx, BAR_Y + BAR_H + 10), (tx, BAR_Y + BAR_H + 32)],
                      fill=(*AMBER, la), width=2)

    # Introductory annotation
    if t > 0.55:
        ann_a = int(220 * (t - 0.55) / 0.45)
        ann   = "NNK degenerate codons introduced at target sites"
        tw    = _textw(F_MED, ann)
        draw.text((WIDTH//2 - tw//2, 108), ann, font=F_MED, fill=(*MID, ann_a))

    la = int(255 * _ease_io(max(0.0, (t - 0.58) / 0.42)))
    _scene_text(draw, "Target",
                "Site-saturation mutagenesis diversifies specific codon positions", la)
    return img


# ═══════════════════════════════════════════════════════════════════════════════
# SCENE 2 — DECODE
# ═══════════════════════════════════════════════════════════════════════════════

def scene_decode(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)

    # Static mini reference bar at top
    _draw_roi_bar(img, 110, 38, 980, 22, alpha=200)

    draw = ImageDraw.Draw(img)

    # Column centres for the 3 target sites
    COL_CX   = [250, 600, 950]
    CHIP_W, CHIP_H = 120, 44
    CHIP_SP  = 54   # spacing between chips in a column

    for col_i, (cx, site_obs) in enumerate(zip(COL_CX, TOP6)):
        # Column header
        col_a = int(255 * min(1.0, t * 5))
        hdr   = TARGET_SITE_NAMES[col_i]
        tw    = _textw(F_SM, hdr)
        draw.text((cx - tw//2, 70), hdr, font=F_SM, fill=(*AMBER, col_a))

        for rank, (aa, freq) in enumerate(site_obs):
            col_delay  = col_i * 0.06
            rank_delay = rank  * 0.055
            total_delay = col_delay + rank_delay
            local_t    = max(0.0, min(1.0,
                            (t - total_delay - 0.08) / max(1e-4, 1.0 - total_delay - 0.08)))
            anim_t     = _elastic_out(local_t)
            a          = int(255 * min(1.0, local_t * 4))
            if a < 5:
                continue

            y_final = 92 + rank * CHIP_SP
            y       = int(y_final - 60 * (1.0 - anim_t))
            chip_x  = cx - CHIP_W // 2

            col = _AA_COL.get(aa, SLATE)
            draw.rounded_rectangle((chip_x, y, chip_x+CHIP_W, y+CHIP_H),
                                   radius=10, fill=(*col, a))

            # AA letter
            lbl_aa = aa
            tw_aa  = _textw(F_BIG, lbl_aa)
            draw.text((chip_x + CHIP_W//2 - tw_aa//2, y + CHIP_H//2 - 14),
                      lbl_aa, font=F_BIG, fill=(255,255,255,a))

            # Frequency %
            if local_t > 0.5:
                fa  = int(a * (local_t - 0.5) / 0.5)
                pct = f"{freq*100:.1f}%"
                tw_pct = _textw(F_XS, pct)
                draw.text((chip_x + CHIP_W//2 - tw_pct//2, y + CHIP_H - 14),
                          pct, font=F_XS, fill=(220,220,220,fa))

    la = int(255 * _ease_io(max(0.0, (t - 0.60) / 0.40)))
    _scene_text(draw, "Decode",
                "Each read is decoded at the target codons to identify the amino acid", la)
    return img


# ═══════════════════════════════════════════════════════════════════════════════
# SCENE 3 — HEATMAP
# ═══════════════════════════════════════════════════════════════════════════════

def scene_heatmap(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)
    draw = ImageDraw.Draw(img)

    N_AA     = len(AA_ORDER)     # 20
    N_SITES  = len(TARGET_SITES) # 3
    CELL_W   = 92
    CELL_H   = 17
    GAP_X    = 24                # gap between site columns
    LABEL_W  = 26                # width for AA labels on left

    total_w  = N_SITES * CELL_W + (N_SITES - 1) * GAP_X + LABEL_W
    HM_X     = (WIDTH - total_w) // 2
    HM_Y     = 68

    ax_a = int(255 * min(1.0, t * 5))

    # Site column headers
    for si, name in enumerate(TARGET_SITE_NAMES):
        cx = HM_X + LABEL_W + si * (CELL_W + GAP_X) + CELL_W // 2
        tw = _textw(F_SM, name)
        draw.text((cx - tw//2, HM_Y - 22), name, font=F_SM, fill=(*AMBER, ax_a))

    # NNK expected indicator label
    if t > 0.85:
        ea  = int(220 * (t - 0.85) / 0.15)
        elbl = "NNK expected →"
        draw.text((HM_X + total_w + 22, HM_Y - 22), elbl, font=F_SM, fill=(*SLATE, ea))

    max_freq = max(
        max(obs.values()) for obs in SSM_OBSERVED
    )

    for ai, aa in enumerate(AA_ORDER):
        # Row reveal: row ai appears when t ≥ ai/N_AA * 0.78
        row_thresh = ai / N_AA * 0.78 + 0.04
        row_t      = _ease_out(max(0.0, (t - row_thresh) / max(1e-4, 0.22)))
        row_a      = int(255 * row_t)
        if row_a < 4:
            continue

        ry = HM_Y + ai * CELL_H

        # AA label
        draw.text((HM_X, ry + 2), aa, font=F_XS, fill=(*INK, row_a))

        # Expected NNK tick mark (rightmost column side)
        exp_frac = NNK_EXPECTED[aa] / max_freq
        exp_x    = HM_X + LABEL_W + N_SITES * (CELL_W + GAP_X) - GAP_X + 18
        exp_bar  = max(1, int(14 * exp_frac))
        if t > 0.85:
            ea = int(row_a * (t - 0.85) / 0.15)
            draw.rectangle((exp_x, ry+3, exp_x+exp_bar, ry+CELL_H-3),
                           fill=(*SLATE, ea))

        # Cells for each site
        for si, obs in enumerate(SSM_OBSERVED):
            cx = HM_X + LABEL_W + si * (CELL_W + GAP_X)
            freq  = obs[aa]
            color = _heatmap_color(freq, max_freq)
            draw.rectangle((cx, ry, cx+CELL_W-2, ry+CELL_H-1),
                           fill=(*color, row_a))

            # Frequency label inside cell if significant
            if freq > 0.035 and CELL_W > 30 and row_a > 80:
                lbl = f"{freq*100:.0f}%"
                tw  = _textw(F_XS, lbl)
                tx  = cx + CELL_W//2 - tw//2 - 1
                fa  = min(row_a, int(210 * min(1.0, (freq - 0.035) / 0.06)))
                draw.text((tx, ry + 3), lbl, font=F_XS, fill=(*INK, fa))

    la = int(255 * _ease_io(max(0.0, (t - 0.62) / 0.38)))
    _scene_text(draw, "Heatmap",
                "Observed amino-acid frequency at every target site", la)
    return img


# ═══════════════════════════════════════════════════════════════════════════════
# SCENE 4 — ASSESS
# ═══════════════════════════════════════════════════════════════════════════════

def scene_assess(fi: int) -> Image.Image:
    t   = _ease_io(fi / (SCENE_FRAMES - 1))
    img = _base_img()
    _draw_card(img, 65, 30, 1135, 472)
    draw = ImageDraw.Draw(img)

    # Chart bounds
    CX, CB, CW, CH = 220, 420, 740, 290
    CY    = CB - CH   # = 130
    MAX_V = max(LOAD_VALS) * 1.18
    N     = len(LOAD_VALS)
    SLOT  = CW // N         # px per bar slot
    BAR_W = SLOT - 30

    ax_a = int(255 * min(1.0, t * 5))
    draw.line([(CX, CY-4), (CX, CB+2)], fill=(*MID, ax_a), width=2)
    draw.line([(CX-2, CB), (CX+CW, CB)], fill=(*MID, ax_a), width=2)

    # Y ticks
    for frac, lbl in [(0.1,"10 %"), (0.2,"20 %"), (0.3,"30 %"), (0.4,"40 %")]:
        ty = int(CB - frac / MAX_V * CH)
        draw.line([(CX-5, ty), (CX+3, ty)], fill=(*MID, ax_a), width=1)
        tw = _textw(F_SM, lbl)
        draw.text((CX - tw - 7, ty - 8), lbl, font=F_SM, fill=(*MID, ax_a))

    draw.text((CX-10, CY-22), "fraction of reads", font=F_SM, fill=(*MID, ax_a))
    xlab = "target sites mutated per read"
    draw.text((CX + CW//2 - _textw(F_SM,xlab)//2, CB+10),
              xlab, font=F_SM, fill=(*MID, ax_a))

    for k, (val, lbl, col) in enumerate(zip(LOAD_VALS, LOAD_LABELS, LOAD_COLS)):
        delay = 0.05 * k
        bar_t = _ease_out(max(0.0, (t - delay - 0.05) / 0.95))
        bx    = CX + k * SLOT + 15
        bar_h = int(CH * val / MAX_V * bar_t)

        if bar_h > 0:
            # Highlight the "fully loaded" bar
            is_full = (k == len(LOAD_VALS) - 1)
            outline_col = GREEN if is_full else col
            draw.rounded_rectangle((bx, CB-bar_h, bx+BAR_W, CB),
                                   radius=6, fill=(*col, 210),
                                   outline=(*outline_col, 255 if is_full else 0),
                                   width=3 if is_full else 0)

        if bar_t > 0.6:
            ta      = int(220 * (bar_t - 0.6) / 0.4)
            pct_str = f"{val*100:.0f}%"
            tw      = _textw(F_SM, pct_str)
            draw.text((bx + BAR_W//2 - tw//2, CB-bar_h-22),
                      pct_str, font=F_SM, fill=(*INK, ta))

        draw.text((bx + BAR_W//2 - _textw(F_SM,lbl)//2, CB+8),
                  lbl, font=F_SM, fill=(*MID, ax_a))

    # Summary annotation
    if t > 0.65:
        ta  = int(255 * (t - 0.65) / 0.35)
        ann = "92% of reads carry ≥ 1 mutation  ·  mean load = 1.8 sites / read"
        tw  = _textw(F_MED, ann)
        draw.text((WIDTH//2 - tw//2, CY + 8), ann, font=F_MED, fill=(*BLUE, ta))

    la = int(255 * _ease_io(max(0.0, (t - 0.58) / 0.42)))
    _scene_text(draw, "Assess",
                "Mutational load confirms efficient diversification across all target sites", la)
    return img


# ═══════════════════════════════════════════════════════════════════════════════
# ASSEMBLE
# ═══════════════════════════════════════════════════════════════════════════════

def _crossfade(a: Image.Image, b: Image.Image, t: float) -> Image.Image:
    return Image.blend(a.convert("RGBA"), b.convert("RGBA"), _ease_io(t))


def build_frames() -> List[Image.Image]:
    fns      = [scene_target, scene_decode, scene_heatmap, scene_assess]
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


_STATIC_COPY = Path("src/uht_tooling/web/static/animations/ssm_profiler.gif")


def main() -> None:
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    frames = build_frames()
    scene_assess(SCENE_FRAMES - 1).convert("RGB").save(POSTER_PATH)
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
