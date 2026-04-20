#!/usr/bin/env python3
"""Generate the JL monogram favicon set from a single source design.

Outputs (into assets/):
    favicon-32.png        # standard browser tab
    favicon-180.png       # apple-touch-icon
    favicon-512.png       # PWA / large preview
    favicon.ico           # legacy (16, 32, 48 packed)
    og-image.png          # Open Graph card (1200x630)

Run:  .venv/bin/python3 scripts/make_favicon.py
"""

from PIL import Image, ImageDraw, ImageFont
from pathlib import Path

ASSETS = Path(__file__).resolve().parent.parent / "assets"
ASSETS.mkdir(exist_ok=True)

# Site palette
BG = (253, 253, 251)        # --color-bg #fdfdfb
INK = (139, 26, 26)         # --color-accent #8b1a1a
MUTED = (85, 85, 85)        # --color-text-muted #555
DARK = (26, 26, 26)         # --color-text #1a1a1a

FONT_BOLD = "/System/Library/Fonts/Supplemental/Georgia Bold.ttf"
FONT_REG = "/System/Library/Fonts/Supplemental/Georgia.ttf"
FONT_CJK = "/System/Library/Fonts/Supplemental/Songti.ttc"  # has SC + bold cuts


def _fit_font_size(text: str, font_path: str, target_w: int, target_h: int,
                   start: int) -> ImageFont.FreeTypeFont:
    """Pick the largest font size that fits text inside (target_w, target_h)."""
    size = start
    while size > 6:
        f = ImageFont.truetype(font_path, size)
        bbox = f.getbbox(text)
        w, h = bbox[2] - bbox[0], bbox[3] - bbox[1]
        if w <= target_w and h <= target_h:
            return f
        size -= 1
    return ImageFont.truetype(font_path, 6)


def make_monogram(size: int, radius_frac: float = 0.16) -> Image.Image:
    """Render a square monogram tile with rounded corners."""
    img = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)

    radius = int(size * radius_frac)
    draw.rounded_rectangle([(0, 0), (size - 1, size - 1)],
                           radius=radius, fill=BG + (255,))

    # Reserve breathing room on all sides (the round-trip through Pillow's
    # textbbox is ink-tight, so we leave ~12% margin so the glyphs don't kiss
    # the corners or the rounded edge).
    text = "JL"
    margin = max(2, int(size * 0.14))
    target_w = size - 2 * margin
    target_h = size - 2 * margin
    font = _fit_font_size(text, FONT_BOLD, target_w, target_h,
                          start=int(size * 0.85))

    bbox = draw.textbbox((0, 0), text, font=font)
    tw = bbox[2] - bbox[0]
    th = bbox[3] - bbox[1]
    # textbbox returns the ink box; offset accounts for the bbox origin.
    x = (size - tw) // 2 - bbox[0]
    y = (size - th) // 2 - bbox[1]

    # Optical baseline nudge: serif caps look bottom-heavy when geometrically centered.
    y -= int(size * 0.015)

    draw.text((x, y), text, fill=INK + (255,), font=font)
    return img


def save_png(img: Image.Image, name: str) -> None:
    out = ASSETS / name
    img.save(out, "PNG", optimize=True)
    print(f"wrote {out.relative_to(ASSETS.parent)} ({img.size[0]}x{img.size[1]})")


def make_ico() -> None:
    """Multi-resolution .ico for legacy browsers and pinned tabs.

    Pillow's `append_images` for ICO is finicky across versions -- the safe
    pattern is to render one large source image and pass the desired sizes
    via `sizes=`, which Pillow then downscales internally. We render the
    source at 64x64 (the size browsers most often request as a high-DPI
    fallback) so the smaller sizes are derived from a clean original.
    """
    sizes = [(16, 16), (32, 32), (48, 48), (64, 64)]
    src = make_monogram(64, radius_frac=0.16)
    out = ASSETS / "favicon.ico"
    src.save(out, format="ICO", sizes=sizes)
    print(f"wrote {out.relative_to(ASSETS.parent)} ({'+'.join(str(s[0]) for s in sizes)})")


def make_og_card() -> None:
    """1200x630 Open Graph / Twitter card.

    Layout: monogram tile on the left, name + tagline on the right.
    """
    W, H = 1200, 630
    img = Image.new("RGB", (W, H), BG)
    draw = ImageDraw.Draw(img)

    # Thin top accent bar in academic red.
    draw.rectangle([(0, 0), (W, 8)], fill=INK)

    # Monogram tile on the left.
    tile_size = 280
    tile = make_monogram(tile_size, radius_frac=0.10)
    # Re-render with a thin border so it reads as a "tile" against the cream bg.
    bordered = Image.new("RGBA", (tile_size, tile_size), (0, 0, 0, 0))
    bd = ImageDraw.Draw(bordered)
    bd.rounded_rectangle([(0, 0), (tile_size - 1, tile_size - 1)],
                         radius=int(tile_size * 0.10),
                         fill=BG + (255,),
                         outline=(228, 225, 216, 255), width=2)
    bordered.alpha_composite(tile)
    img.paste(bordered, (110, (H - tile_size) // 2), bordered)

    # Right column text.
    text_x = 460
    name_font = ImageFont.truetype(FONT_BOLD, 86)
    cjk_font = ImageFont.truetype(FONT_CJK, 44)
    role_font = ImageFont.truetype(FONT_REG, 34)
    url_font = ImageFont.truetype(FONT_REG, 28)

    draw.text((text_x, 195), "Jihao Liu", fill=DARK, font=name_font)
    draw.text((text_x, 305), "刘济豪", fill=MUTED, font=cjk_font)
    draw.text((text_x, 385),
              "Algebraic geometry · Peking University",
              fill=MUTED, font=role_font)
    draw.text((text_x, H - 130), "jihaoliu.org", fill=INK, font=url_font)

    out = ASSETS / "og-image.png"
    img.save(out, "PNG", optimize=True)
    print(f"wrote {out.relative_to(ASSETS.parent)} ({W}x{H})")


def main() -> None:
    save_png(make_monogram(32), "favicon-32.png")
    save_png(make_monogram(180, radius_frac=0.20), "favicon-180.png")
    save_png(make_monogram(512, radius_frac=0.18), "favicon-512.png")
    make_ico()
    make_og_card()


if __name__ == "__main__":
    main()
