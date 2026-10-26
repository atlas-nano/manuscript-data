"""Resolves the data, output and optional panel-image locations used by every figure script."""
import os
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
EXP = DATA / "experimental"
SIM = DATA / "simulation"
TABLES = DATA / "tables"
OUT = pathlib.Path(os.environ.get("FIGURE_OUTDIR", str(ROOT / "figures" / "out")))
OUT.mkdir(parents=True, exist_ok=True)
PANELS = pathlib.Path(os.environ["FIGURE_PANEL_IMAGES"]) if os.environ.get("FIGURE_PANEL_IMAGES") else None


def data(*parts):
    p = DATA.joinpath(*parts)
    if p.exists():
        return p
    g = p.with_name(p.name + ".gz")
    if g.exists():
        return g
    raise FileNotFoundError(p)


def panel_image(name):
    if PANELS is None:
        return None
    p = PANELS / name
    return p if p.exists() else None


def register_fonts():
    from matplotlib import font_manager
    for d in os.environ.get("FIGURE_FONTS", "").split(os.pathsep):
        if d and os.path.isdir(d):
            for f in ("arial.ttf", "arialbd.ttf", "ariali.ttf", "arialbi.ttf"):
                p = os.path.join(d, f)
                if os.path.exists(p):
                    font_manager.fontManager.addfont(p)
