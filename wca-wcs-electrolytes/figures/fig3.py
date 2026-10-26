"""Composes Figure 3 from the panel A to E plate and the panel F plate written by fig3_panels_ae.py and fig3f.py."""
import os
import matplotlib
from matplotlib import font_manager
from PIL import Image, ImageDraw, ImageFont
from paths import OUT

BASE = OUT / "fig3_panels_ae.png"
F_NEW = OUT / "fig3f_plate.png"
im = Image.open(BASE).convert("RGB")
dr = ImageDraw.Draw(im)
dr.rectangle((1440, 858, im.size[0], im.size[1]), fill="white")
f = Image.open(F_NEW).convert("RGB")
fc = f.crop((0, 0, 698, f.size[1]))
im.paste(fc, (1462, 866))
font_file = os.environ.get("FIGURE_LETTER_FONT") or os.path.join(matplotlib.get_data_path(), "fonts", "ttf",
                                                                  "DejaVuSans-Bold.ttf")
font = ImageFont.truetype(font_file, 40)
BOX = {"A": (21, 61, 42, 71), "B": (771, 812, 42, 71), "C": (1461, 1499, 42, 71),
       "D": (21, 62, 810, 839), "E": (771, 811, 810, 839), "F": (1461, 1492, 810, 839)}
for L, (x0, x1, y0, y1) in BOX.items():
    dr.rectangle((x0 - 3, y0 - 4, x1 + 3, y1 + 4), fill="white")
    dr.text((x0, y1 + 1), L, font=font, fill="black", anchor="ls")
im.save(OUT / "Fig3.png", dpi=(300, 300))
print(f"Fig3: {im.size[0]}x{im.size[1]} px")
