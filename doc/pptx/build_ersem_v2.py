#!/usr/bin/env python
"""
ERSEM schematic builder v2 -- renders ersem_spec.py (the verified logic layer)
into an editable PowerPoint with native shapes and clean orthogonal connectors.

Run:
    ~/mambaforge/envs/xfvcom/bin/python build_ersem_v2.py
Output:
    ERSEM_original_editable.pptx
"""
from pptx import Presentation
from pptx.util import Cm, Pt
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE, MSO_CONNECTOR
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.oxml.ns import qn

import ersem_spec as S

SLIDE_W, SLIDE_H = Cm(33.0), Cm(30.0)

COL = {
    "red":   RGBColor(0xCC, 0x33, 0x33),
    "blue":  RGBColor(0x2E, 0x6D, 0xB4),
    "green": RGBColor(0x2E, 0x8B, 0x2E),
    "black": RGBColor(0x20, 0x20, 0x20),
    "grey":  RGBColor(0xA8, 0xA8, 0xA8),
    "text":  RGBColor(0x10, 0x10, 0x10),
    "white": RGBColor(0xFF, 0xFF, 0xFF),
}
# node-kind styling: (fill, line, text-colour, bold, oval?)
KIND = {
    "gas":     (RGBColor(0xFF, 0xFF, 0xFF), COL["blue"],  COL["text"],  False, True),
    "spec":    (RGBColor(0x1C, 0x1C, 0x1C), RGBColor(0x1C, 0x1C, 0x1C), COL["white"], False, False),
    "capsule": (RGBColor(0xF2, 0xF2, 0xF2), COL["grey"],  COL["text"],  False, True),
    "phyto":   (RGBColor(0x7C, 0xA9, 0x3A), RGBColor(0x4F, 0x6E, 0x22), COL["white"], True, True),
    "bact":    (RGBColor(0xD8, 0x8C, 0x3A), RGBColor(0xA9, 0x66, 0x1E), COL["white"], True, True),
    "zoo":     (RGBColor(0x7E, 0x9B, 0xC4), RGBColor(0x53, 0x6E, 0x96), COL["white"], True, True),
    "hetero":  (RGBColor(0x6F, 0xA8, 0xA8), RGBColor(0x47, 0x7C, 0x7C), COL["white"], True, True),
    "pool":    (RGBColor(0xE8, 0xE8, 0xE8), COL["grey"],  COL["text"],  True,  True),
}
BAND_FILL = {
    "atmos":   RGBColor(0xE9, 0xEE, 0xF2),
    "pelagic": RGBColor(0xBE, 0xD2, 0xE0),
    "benthic": RGBColor(0xCF, 0xC2, 0xA0),
}
GROUP_FILL = {
    "Phytoplankton": (RGBColor(0xCF, 0xDD, 0xB0), COL["green"]),
    "Microbes_pel":  (RGBColor(0xDC, 0xDC, 0xDC), COL["grey"]),
    "Microbes_ben":  (RGBColor(0xDC, 0xCF, 0xB0), COL["grey"]),
    "Consumers":     (RGBColor(0xC8, 0xD2, 0xDC), COL["grey"]),
    "Zoobenthos":    (RGBColor(0xC8, 0xD2, 0xDC), COL["grey"]),
}

prs = Presentation()
prs.slide_width, prs.slide_height = SLIDE_W, SLIDE_H
slide = prs.slides.add_slide(prs.slide_layouts[6])
shp = slide.shapes


# ---------------------------------------------------------------------------
def _txt(tf, text, size, bold, color, anchor=MSO_ANCHOR.MIDDLE):
    tf.vertical_anchor = anchor
    tf.word_wrap = True
    for m in ("margin_left", "margin_right", "margin_top", "margin_bottom"):
        setattr(tf, m, Cm(0.03))
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.CENTER
    for i, line in enumerate(text.split("\n")):
        para = p if i == 0 else tf.add_paragraph()
        para.alignment = PP_ALIGN.CENTER
        r = para.add_run()
        r.text = line
        r.font.size = Pt(size)
        r.font.bold = bold
        r.font.color.rgb = color


def add_box(x, y, w, h, *, shape, fill, line, lw, text="", size=10,
            bold=False, tcolor=COL["text"]):
    sp = shp.add_shape(shape, Cm(x), Cm(y), Cm(w), Cm(h))
    if fill is None:
        sp.fill.background()
    else:
        sp.fill.solid(); sp.fill.fore_color.rgb = fill
    if line is None:
        sp.line.fill.background()
    else:
        sp.line.color.rgb = line; sp.line.width = Pt(lw)
    sp.shadow.inherit = False
    if text:
        _txt(sp.text_frame, text, size, bold, tcolor)
    return sp


def add_label(x, y, w, h, text, *, size, bold=True, color=COL["text"],
              align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP):
    tb = shp.add_textbox(Cm(x), Cm(y), Cm(w), Cm(h))
    tf = tb.text_frame; tf.word_wrap = True; tf.vertical_anchor = anchor
    p = tf.paragraphs[0]; p.alignment = align
    r = p.add_run(); r.text = text
    r.font.size = Pt(size); r.font.bold = bold; r.font.color.rgb = color
    return tb


def _arrowhead(conn, end, color):
    ln = conn.line._get_or_add_ln()
    he = ln.makeelement(qn(f"a:{end}"),
                        {"type": "triangle", "w": "med", "len": "med"})
    ln.append(he)


def _dash(conn):
    ln = conn.line._get_or_add_ln()
    d = ln.makeelement(qn("a:prstDash"), {"val": "dash"})
    ln.append(d)


def seg(x1, y1, x2, y2, color, w, *, head=False, tail=False, dash=False):
    c = shp.add_connector(MSO_CONNECTOR.STRAIGHT, Cm(x1), Cm(y1), Cm(x2), Cm(y2))
    c.line.color.rgb = COL[color]; c.line.width = Pt(w)
    c.shadow.inherit = False
    if dash:
        _dash(c)
    if head:
        _arrowhead(c, "tailEnd", color)
    if tail:
        _arrowhead(c, "headEnd", color)
    return c


# ---------------------------------------------------------------------------
# Node geometry helpers -- ports
# ---------------------------------------------------------------------------
def rect(n):
    return (n["cx"] - n["w"] / 2, n["cy"] - n["h"] / 2, n["w"], n["h"])


def port(n, side):
    x, y, w, h = rect(n)
    return {
        "l": (x, y + h / 2), "r": (x + w, y + h / 2),
        "t": (x + w / 2, y), "b": (x + w / 2, y + h),
        "c": (n["cx"], n["cy"]),
    }[side]


def choose_ports(a, b, route):
    """Pick attachment sides based on relative position + requested route."""
    dx, dy = b["cx"] - a["cx"], b["cy"] - a["cy"]
    if route == "h" or (route == "s" and abs(dx) >= abs(dy)):
        sp = "r" if dx >= 0 else "l"
        dp = "l" if dx >= 0 else "r"
    else:
        sp = "b" if dy >= 0 else "t"
        dp = "t" if dy >= 0 else "b"
    return sp, dp


def draw_edge(e):
    # raw coordinate edge (e.g. interface exchange not tied to nodes)
    if e.get("raw"):
        w = 4.5 if e.get("thick") else 1.6
        seg(e["x1"], e["y1"], e["x2"], e["y2"], e["color"], w,
            head=True, tail=e.get("double", False), dash=e.get("dash", False))
        return
    a, b = S.NODES[e["src"]], S.NODES[e["dst"]]
    color = e["color"]
    w = 4.5 if e.get("thick") else 1.6
    route = e.get("route", "s")
    dash = e.get("dash", False)
    double = e.get("double", False)
    sp = e.get("src_port"); dp = e.get("dst_port")
    if not sp or not dp:
        sp, dp = choose_ports(a, b, route)
    x1, y1 = port(a, sp)
    x2, y2 = port(b, dp)

    pts = [(x1, y1)]
    if route == "s" or sp in ("l", "r") and dp in ("l", "r") and abs(y1 - y2) < 0.05:
        pass  # straight
    elif route == "h" or (sp in ("l", "r")):
        midx = (x1 + x2) / 2
        pts += [(midx, y1), (midx, y2)]
    elif route == "v" or (sp in ("t", "b")):
        midy = (y1 + y2) / 2
        pts += [(x1, midy), (x2, midy)]
    pts.append((x2, y2))

    n = len(pts) - 1
    for i in range(n):
        (px1, py1), (px2, py2) = pts[i], pts[i + 1]
        head = (i == n - 1)
        tail = (i == 0 and double)
        seg(px1, py1, px2, py2, color, w, head=head, tail=tail, dash=dash)


# ===========================================================================
# RENDER
# ===========================================================================
M = 0.6
W = 33.0 - 2 * M

# 1. Bands
for key, bd in S.BANDS.items():
    add_box(M, bd["y0"], W, bd["y1"] - bd["y0"], shape=MSO_SHAPE.RECTANGLE,
            fill=BAND_FILL[key], line=COL["grey"], lw=0.75)
add_label(M + 0.3, 0.7, 8, 1.0, "Atmosphere", size=15)
add_label(33.0 - M - 5.2, S.BANDS["pelagic"]["y0"] + 0.15, 4.9, 1.0, "Pelagic",
          size=18, align=PP_ALIGN.RIGHT)
add_label(33.0 - M - 5.2, S.BANDS["benthic"]["y1"] - 1.3, 4.9, 1.0, "Benthic",
          size=18, align=PP_ALIGN.RIGHT)

# 2. Headers + sub-headers
for name, xc in S.HEADERS.items():
    add_label(xc - 3, 0.05, 6, 0.9, name, size=18, align=PP_ALIGN.CENTER)
add_label(M + 0.4, S.BANDS["pelagic"]["y0"] + 0.2, 4, 0.9, "Carbonate\nSystem",
          size=12)
add_label(M + 5.0, S.BANDS["pelagic"]["y0"] + 0.2, 4, 0.9, "Nutrients", size=12)

# 3. Dashed Inorganics|Organics separator
seg(S.SEP_X, S.BANDS["pelagic"]["y0"] + 0.1, S.SEP_X, S.BANDS["benthic"]["y1"] - 0.1,
    "grey", 1.0, dash=True)

# 4. Group containers (behind nodes)
for g in S.GROUPS:
    fill, line = GROUP_FILL[g["name"]]
    add_box(g["x"], g["y"], g["w"], g["h"], shape=MSO_SHAPE.ROUNDED_RECTANGLE,
            fill=fill, line=line, lw=1.0)
    lab = g.get("label", g["name"])
    lp = g.get("label_pos", "tl")
    if lp == "tl":
        add_label(g["x"] + 0.3, g["y"] - 0.02, 6, 0.6, lab, size=12)
    elif lp == "tr":
        add_label(g["x"] + g["w"] - 3.3, g["y"] + 0.1, 3.2, 0.6, lab, size=12,
                  align=PP_ALIGN.RIGHT)
    elif lp == "br":
        add_label(g["x"] + g["w"] - 3.3, g["y"] + g["h"] - 0.7, 3.2, 0.6, lab,
                  size=12, align=PP_ALIGN.RIGHT)

# 5. Edges (under nodes so heads tuck beneath shapes cleanly)
for e in S.EDGES:
    draw_edge(e)

# 6. Nodes
for name, n in S.NODES.items():
    fill, line, tcol, bold, oval = KIND[n["kind"]]
    shape = MSO_SHAPE.OVAL if oval else MSO_SHAPE.ROUNDED_RECTANGLE
    x, y, w, h = rect(n)
    add_box(x, y, w, h, shape=shape, fill=fill, line=line, lw=1.0,
            text=n.get("text", ""), size=10 if n["kind"] != "spec" else 9,
            bold=bold, tcolor=tcol)

# 7. Legend
lx, ly = M + 0.3, S.BANDS["benthic"]["y1"] - 3.6
legend = [("red", "Inorganic nutrient exchange"),
          ("blue", "Carbonate system / gas exchange"),
          ("green", "Organic matter flow (DOM / POM)"),
          ("black", "Predation / grazing"),
          ("grey", "Sinking flux (pelagic→benthic)")]
add_box(lx - 0.15, ly - 0.35, 8.8, 3.4, shape=MSO_SHAPE.RECTANGLE,
        fill=COL["white"], line=COL["grey"], lw=0.75)
add_label(lx, ly - 0.25, 4, 0.5, "Legend", size=11)
for i, (c, t) in enumerate(legend):
    yy = ly + 0.5 + i * 0.56
    seg(lx + 0.1, yy, lx + 1.4, yy, c, 2.4, head=True)
    add_label(lx + 1.6, yy - 0.28, 7, 0.5, t, size=9, bold=False)

prs.save("ERSEM_original_editable.pptx")
print("wrote ERSEM_original_editable.pptx  |  nodes:", len(S.NODES),
      " edges:", len(S.EDGES))
