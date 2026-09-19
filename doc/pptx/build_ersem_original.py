#!/usr/bin/env python
"""
Build an editable PowerPoint reproduction of the original ERSEM schematic
(docs/images/ERSEM.png, after Butenschoen et al. 2016, GMD).

Goal: every box and label is a native PowerPoint shape with editable text,
preserving the conceptual structure (Inorganics/Organics x Pelagic/Benthic,
functional groups, and colour-coded flux arrows) rather than reproducing the
exact graphic styling.

Colour code of arrows (matching the original figure):
    red    : inorganic nutrient exchange (uptake / remineralisation)
    blue   : carbonate system / DIC / gas exchange
    green  : organic matter flow (POM / DOM production and transfer)
    black  : predation / grazing (trophic transfer)
    grey   : sinking flux (pelagic -> benthic coupling)

Run:
    ~/mambaforge/envs/xfvcom/bin/python build_ersem_original.py
Output:
    ERSEM_original_editable.pptx
"""

from pptx import Presentation
from pptx.util import Cm, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE, MSO_CONNECTOR
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.oxml.ns import qn

# ---------------------------------------------------------------------------
# Canvas geometry (cm).  16:9-ish portrait-leaning canvas matching the source.
# ---------------------------------------------------------------------------
SLIDE_W = Cm(33.0)
SLIDE_H = Cm(30.0)

# Palette -------------------------------------------------------------------
COL = {
    "red":   RGBColor(0xCC, 0x33, 0x33),
    "blue":  RGBColor(0x2E, 0x6D, 0xB4),
    "green": RGBColor(0x2E, 0x8B, 0x2E),
    "black": RGBColor(0x20, 0x20, 0x20),
    "grey":  RGBColor(0x9A, 0x9A, 0x9A),

    "phyto":   RGBColor(0x7C, 0xA9, 0x3A),   # green ellipses
    "phyto_bg": RGBColor(0xCF, 0xDD, 0xB0),
    "bact":    RGBColor(0xD8, 0x8C, 0x3A),   # orange ellipses
    "bact_bg": RGBColor(0xD9, 0xD9, 0xD9),
    "zoo":     RGBColor(0x7E, 0x9B, 0xC4),   # blue-grey ellipses (consumers)
    "zoo_teal": RGBColor(0x6F, 0xA8, 0xA8),  # teal (heterotrophs)
    "inorg":   RGBColor(0xF2, 0xF2, 0xF2),   # white-ish nutrient capsules
    "dom":     RGBColor(0xE8, 0xE8, 0xE8),
    "dark":    RGBColor(0x1C, 0x1C, 0x1C),   # carbonate speciation box

    "pel_bg":  RGBColor(0xBE, 0xD2, 0xE0),   # pelagic background band
    "ben_bg":  RGBColor(0xCF, 0xC2, 0xA0),   # benthic background band
    "atmos":   RGBColor(0xE9, 0xEE, 0xF2),
    "text":    RGBColor(0x10, 0x10, 0x10),
    "white":   RGBColor(0xFF, 0xFF, 0xFF),
}

prs = Presentation()
prs.slide_width = SLIDE_W
prs.slide_height = SLIDE_H
slide = prs.slides.add_slide(prs.slide_layouts[6])  # blank
shapes = slide.shapes


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _set_text(tf, text, size=11, bold=False, color=COL["text"], align=PP_ALIGN.CENTER):
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.alignment = align
    run = p.add_run()
    run.text = text
    f = run.font
    f.size = Pt(size)
    f.bold = bold
    f.color.rgb = color
    return p


def box(x, y, w, h, text="", *, shape=MSO_SHAPE.ROUNDED_RECTANGLE,
        fill=None, line=None, line_w=1.0, size=11, bold=False,
        tcolor=COL["text"], anchor=MSO_ANCHOR.MIDDLE):
    sp = shapes.add_shape(shape, Cm(x), Cm(y), Cm(w), Cm(h))
    if fill is None:
        sp.fill.background()
    else:
        sp.fill.solid()
        sp.fill.fore_color.rgb = fill
    if line is None:
        sp.line.fill.background()
    else:
        sp.line.color.rgb = line
        sp.line.width = Pt(line_w)
    sp.shadow.inherit = False
    tf = sp.text_frame
    tf.vertical_anchor = anchor
    tf.margin_left = Cm(0.05)
    tf.margin_right = Cm(0.05)
    tf.margin_top = Cm(0.02)
    tf.margin_bottom = Cm(0.02)
    if text:
        _set_text(tf, text, size=size, bold=bold, color=tcolor)
    return sp


def label(x, y, w, h, text, *, size=14, bold=True, color=COL["text"],
          align=PP_ALIGN.LEFT, italic=False):
    tb = shapes.add_textbox(Cm(x), Cm(y), Cm(w), Cm(h))
    tf = tb.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.alignment = align
    run = p.add_run()
    run.text = text
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = color
    return tb


def _line_dash(connector, dash="dash"):
    ln = connector.line._get_or_add_ln()
    d = ln.find(qn("a:prstDash"))
    if d is None:
        d = ln.makeelement(qn("a:prstDash"), {})
        ln.append(d)
    d.set("val", dash)


def arrow(x1, y1, x2, y2, *, color="black", w=1.5, dash=None,
          head=True, tail=False):
    """Straight connector with optional arrowheads. Coords in cm."""
    cn = shapes.add_connector(MSO_CONNECTOR.STRAIGHT,
                              Cm(x1), Cm(y1), Cm(x2), Cm(y2))
    cn.line.color.rgb = COL[color]
    cn.line.width = Pt(w)
    if dash:
        _line_dash(cn, dash)
    ln = cn.line._get_or_add_ln()
    if head:
        he = ln.makeelement(qn("a:tailEnd"),
                            {"type": "triangle", "w": "med", "len": "med"})
        ln.append(he)
    if tail:
        te = ln.makeelement(qn("a:headEnd"),
                            {"type": "triangle", "w": "med", "len": "med"})
        ln.append(te)
    cn.shadow.inherit = False
    return cn


def elbow(pts, *, color="black", w=1.5, head=True):
    """Poly-line through a list of (x,y) cm points, arrowhead on last seg."""
    for i in range(len(pts) - 1):
        x1, y1 = pts[i]
        x2, y2 = pts[i + 1]
        arrow(x1, y1, x2, y2, color=color, w=w,
              head=(head and i == len(pts) - 2))


# ===========================================================================
# 1. BACKGROUND BANDS  (Atmosphere / Pelagic / Benthic ; Inorganics/Organics)
# ===========================================================================
M = 0.6                       # outer margin
W = 33.0 - 2 * M              # usable width
PEL_TOP, PEL_BOT = 2.4, 16.2
BEN_TOP, BEN_BOT = 16.2, 29.4

# Atmosphere strip
box(M, 0.6, W, 1.7, fill=COL["atmos"], line=COL["grey"], line_w=0.75)
label(M + 0.3, 0.7, 10, 1.0, "Atmosphere", size=15, bold=True)

# Pelagic / Benthic bands
box(M, PEL_TOP, W, PEL_BOT - PEL_TOP, fill=COL["pel_bg"], line=COL["grey"], line_w=0.75)
box(M, BEN_TOP, W, BEN_BOT - BEN_TOP, fill=COL["ben_bg"], line=COL["grey"], line_w=0.75)
label(33.0 - M - 5.0, PEL_TOP + 0.2, 4.6, 1.0, "Pelagic", size=18, bold=True,
      align=PP_ALIGN.RIGHT)
label(33.0 - M - 5.0, BEN_BOT - 1.4, 4.6, 1.0, "Benthic", size=18, bold=True,
      align=PP_ALIGN.RIGHT)

# Top headers: Inorganics | Organics
label(M + 4.0, 0.05, 8, 0.9, "Inorganics", size=18, bold=True, align=PP_ALIGN.CENTER)
label(M + 16.0, 0.05, 10, 0.9, "Organics", size=18, bold=True, align=PP_ALIGN.CENTER)

# Sub-headers inside inorganics
label(M + 0.4, PEL_TOP + 0.2, 5, 0.8, "Carbonate\nSystem", size=12, bold=True)
label(M + 5.4, PEL_TOP + 0.2, 4, 0.8, "Nutrients", size=12, bold=True)

# Dashed separator between Inorganics and Organics
sep_x = M + 10.6
arrow(sep_x, PEL_TOP + 0.1, sep_x, BEN_BOT - 0.1, color="grey", w=1.0,
      dash="dash", head=False)


# ===========================================================================
# 2. INORGANICS  --  carbonate system + nutrients, mirrored Pelagic/Benthic
# ===========================================================================
def carbonate_column(x, y0, tag):
    """Vertical carbonate-system stack. Returns dict of key shape anchors."""
    s = {}
    cap_w, cap_h = 1.7, 1.5
    # pCO2 (top, only pelagic links to atmosphere)
    s["pco2"]  = box(x, y0,        cap_w, cap_h, "pCO₂", shape=MSO_SHAPE.OVAL,
                     fill=COL["white"], line=COL["blue"], size=10)
    s["dic_spec"] = box(x - 0.15, y0 + 2.0, 2.0, 2.4,
                        "H₂CO₃\n↓\nHCO₃⁻\n↓\nCO₃²⁻",
                        shape=MSO_SHAPE.ROUNDED_RECTANGLE, fill=COL["dark"],
                        line=COL["dark"], size=9, tcolor=COL["white"])
    s["ph"]    = box(x + 0.05, y0 + 4.7, 1.6, 1.3, "pH\nΩ", shape=MSO_SHAPE.OVAL,
                     fill=COL["white"], line=COL["blue"], size=9)
    s["ta"]    = box(x - 1.7, y0 + 4.7, 1.5, 1.3, "TA", shape=MSO_SHAPE.OVAL,
                     fill=COL["white"], line=COL["blue"], size=10)
    s["dic"]   = box(x + 0.6, y0 + 1.2, 1.4, 1.2, "DIC", shape=MSO_SHAPE.OVAL,
                     fill=COL["white"], line=COL["blue"], size=9)
    s["o2"]    = box(x + 1.9, y0 + 1.2, 1.4, 1.2, "O₂", shape=MSO_SHAPE.OVAL,
                     fill=COL["white"], line=COL["blue"], size=10)
    return s


def nutrient_column(x, y0, species):
    """Vertical stack of nutrient capsules. Returns dict name->shape."""
    s = {}
    cap_w, cap_h, gap = 1.6, 1.15, 0.18
    for i, name in enumerate(species):
        yy = y0 + i * (cap_h + gap)
        s[name] = box(x, yy, cap_w, cap_h, name, shape=MSO_SHAPE.OVAL,
                      fill=COL["inorg"], line=COL["grey"], size=10)
    return s


# --- Pelagic inorganics
pel_carb = carbonate_column(M + 1.2, PEL_TOP + 1.3, "pel")
pel_nut  = nutrient_column(M + 5.7, PEL_TOP + 1.3,
                           ["Fe", "Si", "PO₄", "NO₃", "NH₄"])

# --- Benthic inorganics (mirror; adds S)
ben_carb = carbonate_column(M + 1.2, BEN_TOP + 0.8, "ben")
ben_nut  = nutrient_column(M + 5.7, BEN_TOP + 0.9,
                           ["S", "PO₄", "NO₃", "NH₄"])

# Atmosphere <-> pelagic pCO2 / O2 gas exchange (blue, double headed)
arrow(M + 2.05, 1.9, M + 2.05, PEL_TOP + 1.3, color="blue", w=1.5, tail=True)


# ===========================================================================
# 3. ORGANICS -- PELAGIC
# ===========================================================================
ORG_X = M + 11.4                      # left edge of organic compartment

# ---- Phytoplankton group box + 4 species ellipses
phx, phy, phw, phh = ORG_X, PEL_TOP + 0.7, 9.2, 2.6
box(phx, phy, phw, phh, fill=COL["phyto_bg"], line=COL["green"], line_w=1.0)
label(phx + 0.3, phy - 0.05, 6, 0.7, "Phytoplankton", size=12, bold=True)
phyto_names = ["Picoph.", "Nanoph.", "Microph.", "Diatoms"]
phyto = {}
ew, eh = 2.0, 1.5
for i, nm in enumerate(phyto_names):
    ex = phx + 0.35 + i * 2.22
    ey = phy + 0.85 - (0.0 if i % 2 else 0.0)
    phyto[nm] = box(ex, ey, ew, eh, nm, shape=MSO_SHAPE.OVAL,
                    fill=COL["phyto"], line=COL["green"], size=10,
                    tcolor=COL["white"], bold=True)

# ---- Microbes (Bacteria)
mbx, mby, mbw, mbh = ORG_X, PEL_TOP + 4.0, 3.2, 2.0
box(mbx, mby, mbw, mbh, fill=COL["bact_bg"], line=COL["grey"], line_w=1.0)
label(mbx + 0.25, mby - 0.05, 3, 0.6, "Microbes", size=11, bold=True)
pel_bact = box(mbx + 0.45, mby + 0.65, 2.3, 1.15, "Bacteria", shape=MSO_SHAPE.OVAL,
               fill=COL["bact"], line=COL["bact"], size=10, tcolor=COL["white"], bold=True)

# ---- Consumers (Heterotrophs, Microzoopl., Mesozoopl.)
csx, csy, csw, csh = ORG_X, PEL_TOP + 6.6, 9.6, 2.3
box(csx, csy, csw, csh, fill=RGBColor(0xC8, 0xD2, 0xDC), line=COL["grey"], line_w=0.75)
label(csx + csw - 3.2, csy + csh - 0.75, 3.0, 0.6, "Consumers", size=11, bold=True,
      align=PP_ALIGN.RIGHT)
hetero = box(csx + 0.4, csy + 0.45, 2.4, 1.4, "Hetero-\ntrophs", shape=MSO_SHAPE.OVAL,
             fill=COL["zoo_teal"], line=COL["zoo_teal"], size=10, tcolor=COL["white"], bold=True)
microz = box(csx + 3.3, csy + 0.45, 2.7, 1.4, "Microzoopl.", shape=MSO_SHAPE.OVAL,
             fill=COL["zoo"], line=COL["zoo"], size=10, tcolor=COL["white"], bold=True)
mesoz  = box(csx + 6.4, csy + 0.45, 2.7, 1.4, "Mesozoopl.", shape=MSO_SHAPE.OVAL,
             fill=COL["zoo"], line=COL["zoo"], size=10, tcolor=COL["white"], bold=True)

# ---- DOM and Particulates (pelagic), right side
domx = ORG_X + 9.7
pel_dom = box(domx, PEL_TOP + 3.3, 2.6, 1.5, "DOM", shape=MSO_SHAPE.OVAL,
              fill=COL["dom"], line=COL["grey"], size=11, bold=True)
pel_part = box(domx, PEL_TOP + 6.0, 2.6, 1.5, "Particulates", shape=MSO_SHAPE.OVAL,
               fill=COL["dom"], line=COL["grey"], size=10, bold=True)


# ===========================================================================
# 4. ORGANICS -- BENTHIC
# ===========================================================================
# ---- Zoobenthos group
zbx, zby, zbw, zbh = ORG_X, BEN_TOP + 0.8, 6.4, 3.6
box(zbx, zby, zbw, zbh, fill=RGBColor(0xC8, 0xD2, 0xDC), line=COL["grey"], line_w=0.75)
label(zbx + zbw - 3.2, zby + 0.1, 3.0, 0.6, "Zoobenthos", size=11, bold=True,
      align=PP_ALIGN.RIGHT)
suspf = box(zbx + 0.4, zby + 0.7, 2.5, 1.4, "Suspension\nFeeders", shape=MSO_SHAPE.OVAL,
            fill=COL["zoo"], line=COL["zoo"], size=10, tcolor=COL["white"], bold=True)
meiob = box(zbx + 0.4, zby + 2.0, 2.5, 1.4, "Meio-\nbenthos", shape=MSO_SHAPE.OVAL,
            fill=COL["zoo"], line=COL["zoo"], size=10, tcolor=COL["white"], bold=True)
depf  = box(zbx + 3.3, zby + 1.4, 2.5, 1.4, "Deposit\nFeeders", shape=MSO_SHAPE.OVAL,
            fill=COL["zoo"], line=COL["zoo"], size=10, tcolor=COL["white"], bold=True)

# ---- Benthic microbes: aerobic + anaerobic bacteria
bmx, bmy, bmw, bmh = ORG_X, BEN_TOP + 4.8, 3.4, 3.4
box(bmx, bmy, bmw, bmh, fill=COL["bact_bg"], line=COL["grey"], line_w=1.0)
label(bmx + bmw - 2.7, bmy + bmh - 0.7, 2.6, 0.6, "Microbes", size=11, bold=True,
      align=PP_ALIGN.RIGHT)
aerob = box(bmx + 0.45, bmy + 0.5, 2.5, 1.15, "Aerobic\nBacteria", shape=MSO_SHAPE.OVAL,
            fill=COL["bact"], line=COL["bact"], size=9, tcolor=COL["white"], bold=True)
anaer = box(bmx + 0.45, bmy + 1.95, 2.5, 1.15, "Anaerobic\nBacteria", shape=MSO_SHAPE.OVAL,
            fill=COL["bact"], line=COL["bact"], size=9, tcolor=COL["white"], bold=True)

# ---- Benthic DOM + Particulates
ben_part = box(domx, BEN_TOP + 3.4, 2.6, 1.5, "Particulates", shape=MSO_SHAPE.OVAL,
               fill=COL["dom"], line=COL["grey"], size=10, bold=True)
ben_dom  = box(domx, BEN_TOP + 7.0, 2.6, 1.5, "DOM", shape=MSO_SHAPE.OVAL,
               fill=COL["dom"], line=COL["grey"], size=11, bold=True)


# ===========================================================================
# 5. ARROWS / FLUXES
# ===========================================================================
# --- Nutrient uptake by phytoplankton (red): nutrients -> phyto group
ny = PEL_TOP + 1.3
arrow(M + 7.3, ny + 1.5, phx, phy + 1.3, color="red", w=1.6)
# --- Remineralisation: bacteria -> nutrients (red)
arrow(mbx, mby + 1.0, M + 7.3, ny + 3.0, color="red", w=1.6)

# --- Carbonate / DIC link to phyto (blue): DIC -> phyto
arrow(M + 3.0, PEL_TOP + 2.5, phx, phy + 1.8, color="blue", w=1.5)

# --- Phytoplankton -> DOM and Particulates (green)
elbow([(phx + phw, phy + 1.0), (domx + 1.3, phy + 1.0), (domx + 1.3, PEL_TOP + 3.3)],
      color="green")
elbow([(phx + phw, phy + 2.0), (domx + 0.6, phy + 2.0), (domx + 0.6, PEL_TOP + 6.0)],
      color="green")
# --- Bacteria <-> DOM (green, uptake) ; DOM -> bacteria
arrow(pel_dom.left / 360000, PEL_TOP + 4.05, mbx + mbw, mby + 1.0, color="green", w=1.5)

# --- Grazing (black): phyto -> consumers ; microbial loop
arrow(phx + 4.0, phy + phh, microz.left/360000 + 1.3, csy + 0.45, color="black", w=1.4)
arrow(pel_bact.left/360000 + 1.1, mby + mbh, hetero.left/360000 + 1.2, csy + 0.45,
      color="black", w=1.4)
arrow(hetero.left/360000 + 2.4, csy + 1.1, microz.left/360000, csy + 1.1, color="black", w=1.4)
arrow(microz.left/360000 + 2.7, csy + 1.1, mesoz.left/360000, csy + 1.1, color="black", w=1.4)
# consumers -> particulates / DOM (green, egestion/mortality)
elbow([(mesoz.left/360000 + 2.7, csy + 1.1), (domx + 1.3, csy + 1.1),
       (domx + 1.3, PEL_TOP + 7.5)], color="green")

# --- Sinking flux pelagic Particulates -> benthic Particulates (grey, thick)
arrow(domx + 1.3, PEL_TOP + 7.5, domx + 1.3, BEN_TOP + 3.4, color="grey", w=5.0)

# --- Benthic: particulates -> DOM, -> bacteria, bacteria -> nutrients (red/green)
arrow(domx, BEN_TOP + 4.1, bmx + bmw, bmy + 1.0, color="green", w=1.5)
arrow(bmx, bmy + 1.5, M + 7.3, BEN_TOP + 2.0, color="red", w=1.6)
elbow([(ben_part.left/360000, BEN_TOP + 4.9), (domx + 1.3, BEN_TOP + 6.4),
       (domx + 1.3, BEN_TOP + 7.0)], color="green")
# zoobenthos grazing on particulates/bacteria (black)
arrow(zbx + zbw, zby + 1.8, ben_part.left/360000, BEN_TOP + 4.1, color="black", w=1.4)

# Benthic carbonate/O2 mirror gas exchange to pelagic (blue dashed up the column)
arrow(M + 2.05, BEN_TOP + 0.8, M + 2.05, PEL_BOT, color="blue", w=1.2,
      dash="dash", head=False)

# --- Legend -----------------------------------------------------------------
lx, ly = M + 0.3, BEN_BOT - 3.4
legend = [
    ("red",   "Inorganic nutrient exchange"),
    ("blue",  "Carbonate system / gas exchange"),
    ("green", "Organic matter flow (DOM / POM)"),
    ("black", "Predation / grazing"),
    ("grey",  "Sinking flux (pelagic→benthic)"),
]
box(lx - 0.15, ly - 0.3, 8.6, 3.3, fill=COL["white"], line=COL["grey"], line_w=0.75)
label(lx, ly - 0.2, 4, 0.5, "Legend", size=11, bold=True)
for i, (c, txt) in enumerate(legend):
    yy = ly + 0.45 + i * 0.55
    arrow(lx + 0.1, yy + 0.12, lx + 1.3, yy + 0.12, color=c, w=2.2)
    label(lx + 1.5, yy - 0.12, 7, 0.5, txt, size=9, bold=False)

prs.save("ERSEM_original_editable.pptx")
print("wrote ERSEM_original_editable.pptx")
