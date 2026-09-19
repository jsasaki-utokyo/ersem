#!/usr/bin/env python
"""
Render a .pptx to PNG by reading native shapes with python-pptx and drawing
them with matplotlib.  Used to visually verify the generated schematic when
LibreOffice is not available.

Usage:
    ~/mambaforge/envs/xfvcom/bin/python preview_pptx.py <file.pptx> [out.png]
"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse, Rectangle, FancyArrow
from pptx import Presentation
from pptx.util import Emu
from pptx.enum.shapes import MSO_SHAPE_TYPE

EMU_PER_CM = 360000.0


def cm(v):
    return v / EMU_PER_CM


def rgb_of(color_fmt, default=None):
    try:
        c = color_fmt.rgb
        return "#%02x%02x%02x" % (c[0], c[1], c[2])
    except Exception:
        return default


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "ERSEM_original_editable.pptx"
    out = sys.argv[2] if len(sys.argv) > 2 else path.replace(".pptx", "_preview.png")
    prs = Presentation(path)
    sw, sh = cm(prs.slide_width), cm(prs.slide_height)
    slide = prs.slides[0]

    fig, ax = plt.subplots(figsize=(sw / 2.54, sh / 2.54), dpi=130)
    ax.set_xlim(0, sw)
    ax.set_ylim(0, sh)
    ax.invert_yaxis()
    ax.axis("off")

    for sp in slide.shapes:
        st = sp.shape_type
        # Connectors (arrows / lines)
        if st == MSO_SHAPE_TYPE.LINE or sp.__class__.__name__ == "Connector":
            x1, y1 = cm(sp.begin_x), cm(sp.begin_y)
            x2, y2 = cm(sp.end_x), cm(sp.end_y)
            col = rgb_of(sp.line.color, "#333333")
            lw = (sp.line.width or 12700) / 12700.0  # EMU->pt
            ax.plot([x1, x2], [y1, y2], color=col, lw=lw, solid_capstyle="round",
                    zorder=3)
            # arrowheads
            from pptx.oxml.ns import qn
            ln = sp.line._get_or_add_ln()
            if ln.find(qn("a:tailEnd")) is not None:
                ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                            arrowprops=dict(arrowstyle="-|>", color=col,
                                            lw=lw, mutation_scale=12),
                            zorder=6)
            if ln.find(qn("a:headEnd")) is not None:
                ax.annotate("", xy=(x1, y1), xytext=(x2, y2),
                            arrowprops=dict(arrowstyle="-|>", color=col,
                                            lw=lw, mutation_scale=12),
                            zorder=6)
            continue

        x, y, w, h = cm(sp.left), cm(sp.top), cm(sp.width), cm(sp.height)
        fill = None
        try:
            if sp.fill.type is not None and sp.fill.type == 1:  # solid
                fill = rgb_of(sp.fill.fore_color, None)
        except Exception:
            fill = None
        edge = rgb_of(sp.line.color, "none")

        is_oval = False
        try:
            if sp.auto_shape_type is not None and "OVAL" in str(sp.auto_shape_type):
                is_oval = True
        except Exception:
            pass

        if is_oval:
            patch = Ellipse((x + w / 2, y + h / 2), w, h,
                            facecolor=fill or "none", edgecolor=edge, lw=1, zorder=2)
        else:
            patch = FancyBboxPatch((x, y), w, h,
                                   boxstyle="round,pad=0,rounding_size=0.15",
                                   facecolor=fill or "none", edgecolor=edge, lw=1,
                                   zorder=1 if fill else 2)
        ax.add_patch(patch)

        if sp.has_text_frame and sp.text_frame.text.strip():
            txt = sp.text_frame.text
            # crude colour/bold pick from first run
            tcol = "#101010"
            bold = False
            size = 10
            for para in sp.text_frame.paragraphs:
                for run in para.runs:
                    if run.font.color and run.font.color.type is not None:
                        tcol = rgb_of(run.font.color, "#101010")
                    bold = bool(run.font.bold)
                    if run.font.size:
                        size = run.font.size.pt
                    break
                break
            ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center",
                    fontsize=max(5, size * 0.62), color=tcol,
                    fontweight="bold" if bold else "normal", zorder=5,
                    linespacing=0.95)

    fig.tight_layout(pad=0.1)
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
