#!/usr/bin/env python3
"""Kronrod — GIF artifact generator (deterministic).

Rasterizes the 16 committed SVG frames (assets/frames.svg — the adaptive
bisection of the unit quarter disc, (10,21) rule) into:

  assets/anim.gif  — 16-frame animation, 0.5 s per frame, infinite loop
  assets/sheet.gif — static 4x4 contact sheet (0.5 scale, same layout)

The frame geometry (interval bars, bar heights, disc, captions) is parsed
verbatim from the committed SVG, so the GIF is a pixel-faithful raster of
assets/frames.svg. No timestamps, no randomness: the output bytes are a pure
function of assets/frames.svg and the palette below. Re-run any time; the
hashes must match.
"""

import html
import os
import re

from PIL import Image, ImageDraw, ImageFont

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ASSETS = os.path.join(ROOT, "assets")
FRAMES = os.path.join(ASSETS, "frames.svg")

# House palette (src/svg.rs).
PAPER = (0xF7, 0xF6, 0xF1)
INK = (0x25, 0x25, 0x24)
MUTED = (0x67, 0x66, 0x62)
HAIR = (0xDC, 0xDA, 0xD1)
GREEN = (0x22, 0xAC, 0x80)
INDIGO = (0x5B, 0x51, 0xC7)
ORANGE = (0xE8, 0x83, 0x3A)
BLACK = (0, 0, 0)
WHITE = (255, 255, 255)

PALETTE_COLORS = [PAPER, INK, MUTED, HAIR, GREEN, INDIGO, ORANGE, BLACK, WHITE]
FILLS = {
    "#F7F6F1": PAPER,
    "#252524": INK,
    "#676662": MUTED,
    "#DCDAD1": HAIR,
    "#22AC80": GREEN,
    "#5B51C7": INDIGO,
    "#E8833A": ORANGE,
}

W, H = 320, 256


def fonts():
    base = r"C:\Windows\Fonts"
    reg = os.path.join(base, "cour.ttf")
    bold = os.path.join(base, "courbd.ttf")
    cache = {}

    def f(size, b=False):
        key = (int(size), b)
        if key not in cache:
            cache[key] = ImageFont.truetype(bold if b else reg, int(size))
        return cache[key]

    return f


def to_palette(img):
    pal = Image.new("P", (1, 1), 0)
    flat = [c for col in PALETTE_COLORS for c in col]
    flat += [0] * (768 - len(flat))
    pal.putpalette(flat)
    return img.convert("P", palette=pal, dither=Image.Dither.NONE)


def blend(fg, bg, t):
    """fg over bg with opacity t (SVG opacity semantics)."""
    return tuple(int(round(fg[i] * t + bg[i] * (1.0 - t))) for i in range(3))


def cells(svg_text):
    starts = [m.end() for m in re.finditer(r'<g transform="translate\([\d.]+,[\d.]+\)">', svg_text)]
    out = [svg_text[starts[i]:starts[i + 1]] for i in range(len(starts) - 1)]
    out.append(svg_text[starts[-1]:])
    return [ch[:ch.index("</g>")] if "</g>" in ch else ch for ch in out]


def parse_cell(ch):
    """One frame: disc geometry, interval bars, texts."""
    m = re.search(r'<path d="M ([\d.]+) ([\d.]+) L ([\d.]+) ([\d.]+) A ([\d.]+) ([\d.]+)', ch)
    ox, oy, r = float(m.group(1)), float(m.group(2)), float(m.group(5))
    lm = re.search(r'<line x1="([\d.]+)" y1="([\d.]+)" x2="([\d.]+)" y2="([\d.]+)"', ch)
    axis = (float(lm.group(1)), float(lm.group(2)), float(lm.group(3)), float(lm.group(4)))
    bars = [
        (float(a), float(b), float(c), float(d))
        for a, b, c, d in re.findall(r'<rect x="([\d.]+)" y="([\d.]+)" width="([\d.]+)" height="([\d.]+)" fill="#5B51C7"', ch)
    ]
    texts = []
    for x, y, size, fill, content in re.findall(
        r'<text x="([\d.]+)" y="([\d.]+)"[^>]*font-size="([\d.]+)"[^>]*fill="(#[0-9A-Fa-f]{6})"[^>]*>([^<]*)', ch
    ):
        texts.append((float(x), float(y), float(size), FILLS[fill], html.unescape(content)))
    return ox, oy, r, axis, bars, texts


def render_frame(cell, f):
    ox, oy, r, axis, bars, texts = parse_cell(cell)
    img = Image.new("RGB", (W, H), PAPER)
    d = ImageDraw.Draw(img)
    bbox = [ox - r, oy - r, ox + r, oy + r]
    # quarter disc: faint fill + arc (the integrand)
    d.pieslice(bbox, 270, 360, fill=blend(INDIGO, PAPER, 0.07))
    d.arc(bbox, 270, 360, fill=INK, width=2)
    # axis
    d.line([(axis[0], axis[1]), (axis[2], axis[3])], fill=HAIR, width=1)
    # error bars below the axis
    for x, y, w, h in bars:
        d.rectangle([x, y, x + w, y + h], fill=INDIGO)
    # captions
    for x, y, size, color, s in texts:
        font = f(size)
        d.text((x, y - 0.8 * size), s, font=font, fill=color)
    return img


def make_anim(path, frames):
    frames[0].save(
        path,
        save_all=True,
        append_images=frames[1:],
        duration=500,
        loop=0,
        optimize=False,
    )


def make_sheet(path, frames, title):
    scale = 0.5
    cw, ch = int(W * scale), int(H * scale)  # 160 x 128
    gap, margin, title_h = 8, 8, 24
    cols, rows = 4, 4
    w = margin * 2 + gap * (cols - 1) + cols * cw
    h = title_h + margin + gap * (rows - 1) + rows * ch
    img = Image.new("RGB", (w, h), PAPER)
    d = ImageDraw.Draw(img)
    f = fonts()
    d.text((8, 4), title, font=f(11), fill=GREEN)
    for i, fr in enumerate(frames):
        r, col = divmod(i, cols)
        x = margin + gap + col * (cw + gap)
        y = title_h + margin + gap + r * (ch + gap)
        small = fr.resize((cw, ch), Image.Resampling.LANCZOS)
        img.paste(small, (x, y))
        d.rectangle([x, y, x + cw - 1, y + ch - 1], outline=HAIR)
    to_palette(img).save(path, optimize=False)


def main():
    s = open(FRAMES, encoding="utf-8").read()
    chs = cells(s)
    assert len(chs) == 16, "expected 16 committed frames, got %d" % len(chs)
    title_m = re.search(r'<text[^>]*>([^<]*quarter disc[^<]*)</text>', s)
    title = html.unescape(title_m.group(1)) if title_m else "Kronrod — adaptive bisection (pi/4)"
    f = fonts()
    frames = [to_palette(render_frame(c, f)) for c in chs]
    os.makedirs(ASSETS, exist_ok=True)
    anim = os.path.join(ASSETS, "anim.gif")
    sheet = os.path.join(ASSETS, "sheet.gif")
    make_anim(anim, frames)
    make_sheet(sheet, frames, title)
    import hashlib

    for p in (anim, sheet):
        h = hashlib.sha256(open(p, "rb").read()).hexdigest()
        print("wrote %s (%d bytes) sha256 %s" % (os.path.basename(p), os.path.getsize(p), h))


if __name__ == "__main__":
    main()
